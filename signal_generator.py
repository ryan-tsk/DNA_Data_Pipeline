from fast5_research import Fast5
import scrappy
from uuid import uuid4

import numpy as np
import os


def generate_signal(data, file_prefix, directory):

    for i, seq in enumerate(data):
        filename = f'{file_prefix}{i}.fast5'
        out_filepath = os.path.join(directory, filename)
        simulate_read(seq, out_filepath)

    return directory


def simulate_read(seq, out_filepath):
    """
    Simulates nanopore sequencer, writes to file
    Parameters
    ----------
    seq : string
        Nucleotide sequence
    out_dir : string
        Working directory to write intermediate results to.

    out_filepath : string
        Name of FAST5 file which will be written as result of read.
        directory by default.
    id : string
        ID given to read
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
        'read_id': str(uuid4()),
        'scaling_used': 1,
        'median_before': 0,
    }
    tracking_id = {
        'exp_start_time': '1970-01-01T00:00:00Z',
        'run_id': str(uuid4()).replace('-', ''),
        'flow_cell_id': 'FAH00000',
        'device_id': 'TEST123',
        'sample_id': 'TEST123',
    }
    context_tags = {}

    with Fast5.New(out_filepath, 'w', tracking_id=tracking_id, context_tags=context_tags, channel_id=channel_id) as h:
        h.set_raw(raw_data, meta=read_id, read_number=1)
