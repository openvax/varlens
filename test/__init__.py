'''
Utility functions for tests.
'''

import sys
import os
import tempfile
from contextlib import contextmanager

import pandas

if sys.version_info[0] < 3:
    from StringIO import StringIO
else:
    from io import StringIO

def data_path(name):
    '''
    Return the absolute path to a file in the test/data directory.
    The name specified should be relative to test/data.
    '''
    return os.path.join(os.path.dirname(__file__), "data", name)

class Capturing(list):
    def __enter__(self):
        self._stdout = sys.stdout
        sys.stdout = self._stringio = StringIO()
        return self

    def __exit__(self, *args):
        self.extend(self._stringio.getvalue().splitlines())
        sys.stdout = self._stdout

def run_and_parse_csv(function, *args):
    with Capturing() as output:
        function(*args)
    try:
        result = pandas.read_csv(StringIO("\n".join(output)))
    except:
        print("Couldn't parse csv. Function: %s. Args: %s.\nOutput:\n%s"
            % (str(function), str(args), "\n".join(output)))
        raise
    return result

@contextmanager
def temp_file(suffix=".csv"):
    fd = tempfile.NamedTemporaryFile(
        suffix=suffix,
        prefix="test_varlens_",
        delete=False)
    filename = fd.name
    fd.close()
    yield filename
    os.unlink(filename)    

def cols_concat(df, columns, delimiter="-"):
    assert df is not None
    zipped = zip(*[df[c] for c in columns])
    return set([delimiter.join(str(item) for item in row) for row in zipped])