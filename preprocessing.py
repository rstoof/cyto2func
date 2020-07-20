import sys

if len(sys.argv) > 1 and sys.argv[1] == 'with-download':
    import download

import decompose_volume
import normalise

#this plots the results:
#import plot
