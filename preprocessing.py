from urllib.request import urlretrieve
import zipfile
import argparse

from decompose_volume import decompose
from normalise import normalise


description = """
Decorrelates flourescence values and estimated cell volume in
flow cytometry data.
"""

default_url = 'https://data.ncl.ac.uk/ndownloader/articles/12073479/versions/1'
parser = argparse.ArgumentParser(description=description, prog='FlowScatt')
parser.add_argument('--download', dest='url', type=str, nargs='?',
                    const=default_url, default='',
                    help='the URL of the zipped data files to retrieve.')
parser.add_argument('--data-dir', type=str, nargs=1,
                    default='./FCS/', metavar='FILENAME',
                    help='the directory in which to store the data.')
parser.add_argument('-f', '--description-file', type=str, nargs=1,
                    default='file_description.csv', metavar='FILENAME',
                    help='the file which describes the fcs file\'s data.')
parser.add_argument('-c', '--channels', nargs=2, metavar='CHANNEL',
                    default=['FSC_H', 'GFP_H'],
                    help='the channels to use for the analysis.')
parser.add_argument('--min-to-back', action='store_true',
                    default=False,
                    help='''If present, normalisation if performed using
                    the minimum flourescence value, not the
                    autofluorescence plasmid''')
parser.add_argument('-o', '--outfile', type=str, nargs=1, metavar='FILENAME',
                    default='standardised.csv',
                    help='the file path to write the processed data to.')


def download_and_extract(url):
    print("starting download, this may take some time")
    urlretrieve(url, "FCS.zip")
    print("files downloaded,starting unpacking")

    with zipfile.ZipFile("FCS.zip", 'r') as zip_ref:
        zip_ref.extractall(DATA_DIR)


if __name__ == '__main__':
    args = parser.parse_args()

    URL = args.url
    DATA_DIR = args.data_dir
    CHANNELS = args.channels
    DESC_FILE = args.description_file
    MIN_TO_BACK = args.min_to_back
    OUTFILE = args.outfile

    if URL:
        download_and_extract(URL)

    with open('volume_decomposed.csv', 'w', newline='') as outfile:
        decompose(outfile, DESC_FILE, DATA_DIR, CHANNELS)

    with open('volume_decomposed.csv') as infile:
        with open(OUTFILE, 'w', newline='') as outfile:
            normalise(infile, outfile, min_to_back=MIN_TO_BACK)
