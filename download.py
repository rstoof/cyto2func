import urllib

url = 'https://data.ncl.ac.uk/ndownloader/articles/12073479/versions/1'

urllib.urlretrieve(url, "FCS.zip")

import zipfile
with zipfile.ZipFile("FCS.zip", 'r') as zip_ref:
    zip_ref.extractall("/FCS")
