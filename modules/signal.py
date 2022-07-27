from fast5_research import Fast5

import uuid
import scrappy
import numpy as np
import os
import h5py


def generate_signal(data, result_directory: str, file_prefix: str = 'TEST', folder: str = 'signals', encode=False):
    directory = os.path.join(result_directory, folder)

    if not os.path.exists(directory):
        os.mkdir(directory)

    for i, seq in enumerate(data):
        filename = f'{file_prefix}_{i}.fast5'
        read_id = f'{file_prefix}_{i}'
        out_filepath = os.path.join(directory, filename)
        simulate_read(seq, out_filepath, read_id)

        if encode:
            file = h5py.File(out_filepath, 'r+')
            set = file['Raw']['Reads']['Read_1']
            read_id = set.attrs['read_id']
            set.attrs['read_id'] = np.string_(read_id)

    return directory


def simulate_read(seq, out_filepath, read_identity):
    """
    Simulates nanopore sequencer, writes to file
    Parameters
    ----------
    seq : string
        Nucleotide sequence
    out_filepath : string
        Name of FAST5 file which will be written as result of read.
        directory by default.
    read_identity : str
    """

    # Using https://nanoporetech.github.io/fast5_research/examples.html as a reference
    squiggle = scrappy.sequence_to_squiggle(seq, rescale=True).data(as_numpy=True)
    raw_data = np.array([])

    for dwell, mean, stdv in squiggle:
        raw_data = np.append(raw_data, np.random.laplace(mean, stdv/np.sqrt(2), int(round(dwell))))

    start, stop = int(min(raw_data - 1)), int(max(raw_data + 1))
    rng = stop - start
    digitisation = 8192.0
    bins = np.arange(start, stop, rng / digitisation)
    # np.int16 is required, the library will refuse to write anything other
    raw_data = np.digitize(raw_data, bins).astype(np.int16)

    # The following are required meta data
    channel_id = {
        'digitisation': digitisation,
        'offset': 0,
        'range': rng,
        'sampling_rate': 4000,
        'channel_number': 1,
        }
    read_id = {
        'start_time': 0,
        'duration': len(raw_data),
        'read_number': 1,
        'start_mux': 1,
        'read_id': read_identity,
        'scaling_used': 1,
        'median_before': 0,
    }
    tracking_id = {
        'exp_start_time': '1970-01-01T00:00:00Z',
        'run_id': read_identity,
        'flow_cell_id': 'FAH00000',
        'device_id': 'TEST123',
        'sample_id': 'TEST123',
    }
    context_tags = {}

    with Fast5.New(out_filepath, 'w', tracking_id=tracking_id, context_tags=context_tags, channel_id=channel_id) as h:
        h.set_raw(raw_data, meta=read_id, read_number=1)
