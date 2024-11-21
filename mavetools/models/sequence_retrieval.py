import socket
import time
import urllib.error
import urllib.parse
import urllib.request
from Bio import Entrez

from mavetools.models.utils import parseFasta

def is_connected():

    """
    Tests if there is an active internet connection.
    """

    try:
        # connect to the host -- tells us if the host is actually
        # reachable
        socket.create_connection(("1.1.1.1", 53))
        return True
    except OSError:
        pass
    return False


def connection_sleep_cycle():

    """
    A waiting routine for short internet outages.
    """

    while not is_connected():
        print('No connection, sleeping a bit and then try again')
        time.sleep(30)


def getUniprotSequence(uniprot_ac, tries=0):

    """
    Fetches sequence data from Uniprot database.

    Parameters
    ----------

    uniprot_ac
        Uniprot accesion identifier of the protein sequence to fetch.

    tries
        A count parameter for the trial and error loop for short time internet disconnections.

    Returns
    -------

    wildtype_sequence
        Amino acid sequence of the querried protein.
    """

    if uniprot_ac is None:
        print(f'Uniprot Ac is None')
        return None

    if not uniprot_ac[0:3] == 'UPI':
        url = 'https://www.uniprot.org/uniprot/%s.fasta' % uniprot_ac
    else:
        url = 'https://www.uniprot.org/uniparc/%s.fasta' % uniprot_ac
    connection_sleep_cycle()
    try:
        request = urllib.request.Request(url)
        response = urllib.request.urlopen(request, timeout=(tries + 1) * 10)
        page = response.read(9000000).decode('utf-8')
    except:
        if tries < 3:
            return getSequence(uniprot_ac, tries=tries + 1)
        else:
            return None

    lines = page.split("\n")

    wildtype_sequences = []

    for line in lines:
        if line == '':
            continue
        if line[0] == '>':
            continue
        wildtype_sequences.append(line)

    wildtype_sequence = ("".join(wildtype_sequences)).replace(" ", "").replace("\n", "")

    if wildtype_sequence == '':
        return None

    return wildtype_sequence

def get_refseq_sequences(refseqs, seq_type='protein'):

    """
    Fetches sequence data from NCBI refseq database.
    Uses the biopython library.

    Parameters
    ----------

    refseqs
        comma-separated string of refseq identifiers

    seq_type
        type of sequence to be retrieved, either 'nucleotide' or 'protein'

    Returns
    -------

    seq_map
        dictionary of {refseq_identifier:sequence}
    """

    Entrez.email = ''

    ret_type = 'fasta_cds_aa'
    if seq_type == 'protein':
        ret_type = 'fasta'

    net_handle = Entrez.efetch(
        db=seq_type, id=refseqs, rettype=ret_type, retmode="text"
    )
    page = net_handle.read()
    net_handle.close()
    right_split = '_prot'
    left_split = '|'
    if seq_type == 'protein':
        right_split = ' '
        left_split = None
    seq_map = parseFasta(page=page, left_split=left_split, right_split=right_split)

    return seq_map
