#!/usr/bin/env python
#Dan Blankenberg

import sys
import os
import tempfile
import shutil
import optparse
from ftplib import FTP
import tarfile
import zipfile
import gzip
import bz2
import subprocess
try:
    # For Python 3.0 and later
    from urllib.request import urlopen
    from io import BytesIO as StringIO
    from io import UnsupportedOperation
except ImportError:
    # Fall back to Python 2's urllib2
    from urllib2 import urlopen
    from StringIO import StringIO
    UnsupportedOperation = AttributeError
from json import loads, dumps


CHUNK_SIZE = 2**20  # 1mb

DATA_TABLE_NAME = 'tx2gene_table'

def cleanup_before_exit( tmp_dir ):
    if tmp_dir and os.path.exists( tmp_dir ):
        shutil.rmtree( tmp_dir )


def stop_err(msg):
    sys.stderr.write(msg)
    sys.exit(1)


def get_dbkey_dbname_id_name( params, dbkey_description=None ):
#    dbkey = params['param_dict']['dbkey_source']['dbkey']
    #TODO: ensure sequence_id is unique and does not already appear in location file
    sequence_id = params['param_dict']['sequence_id']
    if not sequence_id:
        sequence_id = dbkey #uuid.uuid4() generate and use an uuid instead?

#    if params['param_dict']['dbkey_source']['dbkey_source_selector'] == 'new':
#        dbkey_name = params['param_dict']['dbkey_source']['dbkey_name']
#        if not dbkey_name:
#            dbkey_name = dbkey
#    else:
#        dbkey_name = None
    dbkey = params['param_dict']['dbkey']
    dbkey_name = dbkey_description
    sequence_name = params['param_dict']['sequence_name']
    if not sequence_name:
        sequence_name = dbkey_description
        if not sequence_name:
            sequence_name = dbkey
    return dbkey, dbkey_name, sequence_id, sequence_name


def _get_files_in_ftp_path( ftp, path ):
    path_contents = []
    ftp.retrlines( 'MLSD %s' % ( path ), path_contents.append )
    return [ line.split( ';' )[ -1 ].lstrip() for line in path_contents ]


def _get_stream_readers_for_tar( fh, tmp_dir ):
    fasta_tar = tarfile.open( fileobj=fh, mode='r:*' )
    return [x for x in [fasta_tar.extractfile(member) for member in fasta_tar.getmembers()] if x]


def _get_stream_readers_for_zip( fh, tmp_dir ):
    """
    Unpacks all archived files in a zip file.
    Individual files will be concatenated (in _stream_fasta_to_file)
    """
    fasta_zip = zipfile.ZipFile( fh, 'r' )
    rval = []
    for member in fasta_zip.namelist():
        fasta_zip.extract( member, tmp_dir )
        rval.append( open( os.path.join( tmp_dir, member ), 'rb' ) )
    return rval


def _get_stream_readers_for_gzip( fh, tmp_dir ):
    return [ gzip.GzipFile( fileobj=fh, mode='rb') ]


def _get_stream_readers_for_bz2( fh, tmp_dir ):
    return [ bz2.BZ2File( fh.name, 'rb') ]


def convert_to_tx2gene( rscript_gff_to_tx2gene, fasta_filename, file_type, params ):
    if file_type == 'tx2gene':
        return   #no need to extract tx2gene table
    #print file_type
    #If the file is actually a GFF/GTF file then extract the tx2gene
    gff_temp_filename = tempfile.NamedTemporaryFile().name
    shutil.move(fasta_filename, gff_temp_filename)
    args= ['Rscript']
    args.append(rscript_gff_to_tx2gene)
    args.extend(['-x',gff_temp_filename])
    args.extend(['-o',fasta_filename])
    args.extend(['-t',file_type])
    tmp_stderr = tempfile.NamedTemporaryFile( prefix = "tmp-stderr" )
    return_code = subprocess.call( args=args, shell=False, stderr=tmp_stderr.fileno() )
    #return_code = subprocess.call( args=args, shell=False, stderr=None)
    if return_code:
        tmp_stderr.flush()
        tmp_stderr.seek(0)
        print >> sys.stderr, "Error in process call"
        while True:
            chunk = tmp_stderr.read( CHUNK_SIZE )
            if not chunk:
                break
            sys.stderr.write( chunk )
        sys.exit( return_code )
    tmp_stderr.close()



def _download_file(start, fh):
    tmp = tempfile.NamedTemporaryFile()
    tmp.write(start)
    tmp.write(fh.read())
    tmp.flush()
    tmp.seek(0)
    return tmp


def get_stream_reader(fh, tmp_dir):
    """
    Check if file is compressed and return correct stream reader.
    If file has to be downloaded, do it now.
    """
    magic_dict = {
        b"\x1f\x8b\x08": _get_stream_readers_for_gzip,
        b"\x42\x5a\x68": _get_stream_readers_for_bz2,
        b"\x50\x4b\x03\x04": _get_stream_readers_for_zip,
    }
    start_of_file = fh.read(CHUNK_SIZE)
    try:
        fh.seek(0)
    except UnsupportedOperation:  # This is if fh has been created by urlopen
        fh = _download_file(start_of_file, fh)
    for k,v in magic_dict.items():
        if start_of_file.startswith(k):
            return v(fh, tmp_dir)
    try:  # Check if file is tar file
        if tarfile.open(fileobj=StringIO(start_of_file)):
            return _get_stream_readers_for_tar(fh, tmp_dir)
    except tarfile.ReadError:
        pass
    return fh



def add_fasta_to_table(rscript_gff_to_tx2gene, data_manager_dict, fasta_readers, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, params):
    for data_table_name, data_table_entry in _stream_fasta_to_file(rscript_gff_to_tx2gene, fasta_readers, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, params ):
        if data_table_entry:
            _add_data_table_entry( data_manager_dict, data_table_entry, data_table_name )


def download_from_url(rscript_gff_to_tx2gene, data_manager_dict, params, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, tmp_dir ):
    urls = filter( bool, map( lambda x: x.strip(), params['param_dict']['reference_source']['user_url'].split( '\n' ) ) )
    fasta_readers = [ get_stream_reader(urlopen( url ), tmp_dir) for url in urls ]
    add_fasta_to_table(rscript_gff_to_tx2gene,data_manager_dict, fasta_readers, target_directory, dbkey, dbkey_name, sequence_id,sequence_name, params)


def download_from_history(rscript_gff_to_tx2gene, data_manager_dict, params, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, tmp_dir ):
    #TODO: allow multiple FASTA input files
    input_filename = params['param_dict']['reference_source']['input_fasta']
    if isinstance( input_filename, list ):
        fasta_readers = [ get_stream_reader(open(filename, 'rb'), tmp_dir) for filename in input_filename ]
    else:
        fasta_readers = get_stream_reader(open(input_filename), tmp_dir)
    add_fasta_to_table(rscript_gff_to_tx2gene,data_manager_dict, fasta_readers, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, params)


def copy_from_directory(rscript_gff_to_tx2gene, data_manager_dict, params, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, tmp_dir ):
    input_filename = params['param_dict']['reference_source']['fasta_filename']
    create_symlink = params['param_dict']['reference_source']['create_symlink'] == 'create_symlink'
    if create_symlink:
        data_table_entries = _create_symlink( input_filename, target_directory, dbkey, dbkey_name, sequence_id, sequence_name )
    else:
        if isinstance( input_filename, list ):
            fasta_readers = [ get_stream_reader(open(filename, 'rb'), tmp_dir) for filename in input_filename ]
        else:
            fasta_readers = get_stream_reader(open(input_filename), tmp_dir)
        data_table_entries = _stream_fasta_to_file(rscript_gff_to_tx2gene, fasta_readers, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, params )
    for data_table_name, data_table_entry in data_table_entries:
        if data_table_entry:
            _add_data_table_entry( data_manager_dict, data_table_entry, data_table_name )


def _add_data_table_entry( data_manager_dict, data_table_entry, data_table_name ):
    data_manager_dict['data_tables'] = data_manager_dict.get( 'data_tables', {} )
    data_manager_dict['data_tables'][data_table_name] = data_manager_dict['data_tables'].get( DATA_TABLE_NAME, [] )
    data_manager_dict['data_tables'][data_table_name].append( data_table_entry )
    return data_manager_dict


def _stream_fasta_to_file( rscript_gff_to_tx2gene, fasta_stream, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, params, close_stream=True ):
    fasta_base_filename = "%s_tx2gene.tab" % sequence_id
    fasta_filename = os.path.join( target_directory, fasta_base_filename )
    with open( fasta_filename, 'wb+' ) as fasta_writer:

        if isinstance( fasta_stream, list ) and len( fasta_stream ) == 1:
            fasta_stream = fasta_stream[0]

        if isinstance( fasta_stream, list ):
            last_char = None
            for fh in fasta_stream:
                if last_char not in [ None, '\n', '\r', b'\n', b'\r' ]:
                    fasta_writer.write( b'\n' )
                while True:
                    data = fh.read( CHUNK_SIZE )
                    if data:
                        fasta_writer.write( data )
                        last_char = data[-1]
                    else:
                        break
                if close_stream:
                    fh.close()
        else:
            while True:
                data = fasta_stream.read( CHUNK_SIZE )
                if data:
                    fasta_writer.write( data )
                else:
                    break
            if close_stream:
                fasta_stream.close()

    convert_to_tx2gene( rscript_gff_to_tx2gene,fasta_filename, params['param_dict']['file_type'], params )
    return [ ( DATA_TABLE_NAME, dict( value=sequence_id, dbkey=dbkey, name=sequence_name, path=fasta_base_filename ) ) ]


def compute_fasta_length( fasta_file, out_file, keep_first_word=False ):

    infile = fasta_file
    out = open( out_file, 'w')

    fasta_title = ''
    seq_len = 0

    first_entry = True

    for line in open( infile ):
        line = line.strip()
        if not line or line.startswith( '#' ):
            continue
        if line[0] == '>':
            if first_entry == False:
                if keep_first_word:
                    fasta_title = fasta_title.split()[0]
                out.write( "%s\t%d\n" % ( fasta_title[ 1: ], seq_len ) )
            else:
                first_entry = False
            fasta_title = line
            seq_len = 0
        else:
            seq_len += len(line)

    # last fasta-entry
    if keep_first_word:
        fasta_title = fasta_title.split()[0]
    out.write( "%s\t%d\n" % ( fasta_title[ 1: ], seq_len ) )
    out.close()


def _create_symlink( input_filename, target_directory, dbkey, dbkey_name, sequence_id, sequence_name ):
    fasta_base_filename = "%s.fa" % sequence_id
    fasta_filename = os.path.join( target_directory, fasta_base_filename )
    os.symlink( input_filename, fasta_filename )
    return [  ( DATA_TABLE_NAME, dict( value=sequence_id, dbkey=dbkey, name=sequence_name, path=fasta_base_filename ) ) ]


REFERENCE_SOURCE_TO_DOWNLOAD = dict( url=download_from_url, history=download_from_history, directory=copy_from_directory )


def main():
    #Parse Command Line
    parser = optparse.OptionParser()
    parser.add_option( '-d', '--dbkey_description', dest='dbkey_description', action='store', type="string", default=None, help='dbkey_description' )
    parser.add_option( '-b', '--base_dir', dest='base_dir', action='store', type='string', default=None, help='base_dir')
    parser.add_option( '-t', '--type', dest='file_type', action='store', type='string', default=None, help='file_type')
    (options, args) = parser.parse_args()
    filename = args[0]
    #global DATA_TABLE_NAME
    rscript_gff_to_tx2gene=os.path.join( options.base_dir, 'get_tx2gene_table.R')

    #input_type='gff_gtf'
    #if options.file_type != 'gff_gtf':
    # 	file_type='tx2gene'

    params = loads( open( filename ).read() )
    target_directory = params[ 'output_data' ][0]['extra_files_path']
    os.mkdir( target_directory )
    data_manager_dict = {}

    dbkey, dbkey_name, sequence_id, sequence_name = get_dbkey_dbname_id_name( params, dbkey_description=options.dbkey_description ) 

    if dbkey in [ None, '', '?' ]:
        raise Exception( '"%s" is not a valid dbkey. You must specify a valid dbkey.' % ( dbkey ) )

    # Create a tmp_dir, in case a zip file needs to be uncompressed
    tmp_dir = tempfile.mkdtemp()
    #Fetch the input file
    try:
        REFERENCE_SOURCE_TO_DOWNLOAD[ params['param_dict']['reference_source']['reference_source_selector'] ]( rscript_gff_to_tx2gene, data_manager_dict, params, target_directory, dbkey, dbkey_name, sequence_id, sequence_name, tmp_dir)
    finally:
        cleanup_before_exit(tmp_dir)
    #save info to json file
    open( filename, 'wb' ).write( dumps( data_manager_dict ).encode() )

if __name__ == "__main__":
    main()
