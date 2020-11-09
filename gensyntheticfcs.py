from flowio import create_fcs
from numpy.random import normal
from numpy import array


def writefcs(path, channels):
    """Create and write an FCS file with channels specified by the
    channels dictionary.

    The dictionary should contain PNN labels matched with event data,
    the event data is assumed to be synced across all channels,
    i.e. event 1 will be the first element of each channels event
    data, event 2 the second elements, and so on.

    The channel data should be arraylikes of floats and be of the same
    length.
    """

    pnn_labels = list(channels.keys())
    events = zip(*[channels[label] for label in pnn_labels])
    flattened_events = [data for event in events for data in event]
    with open(path, 'wb') as outfile:
        create_fcs(flattened_events, pnn_labels, outfile)


def genlinearlydependentchannels(µ, σ, n, linear_func):
    channel1data = normal(µ, σ, size=n)
    channel2data = array(list(map(linear_func, channel1data)))
    return {"CH-1": channel1data, "CH-2": channel2data}
