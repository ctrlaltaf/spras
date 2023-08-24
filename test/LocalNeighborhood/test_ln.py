import shutil
import sys
from filecmp import cmp
from pathlib import Path

import pytest

from src.local_neighborhood import LocalNeighborhood
from src.util import compare_files

# TODO consider refactoring to simplify the import
# Modify the path because of the - in the directory
SPRAS_ROOT = Path(__file__).parent.parent.parent.absolute()
sys.path.append(str(Path(SPRAS_ROOT, "docker-wrappers", "LocalNeighborhood")))
TEST_DIR = 'test/LocalNeighborhood/'
OUT_FILE_DEFAULT = TEST_DIR+'output/ln-output.txt'
OUT_FILE_BAD = TEST_DIR+'output/ln-output-bad.txt'

class TestLocalNeighborhood:
    """
    Run the local neighborhood algorithm on the example input files and check the output matches the expected output
    """

    def test_ln(self):
        out_path = Path(OUT_FILE_DEFAULT)
        out_path.unlink(missing_ok=True)
        LocalNeighborhood.run(
            nodes=TEST_DIR+'input/ln-nodes.txt',
            network=TEST_DIR+'input/ln-network.txt',
            output_file=OUT_FILE_DEFAULT
        )
        assert out_path.exists()
        expected_file = TEST_DIR + 'expected_output/ln-output.txt'
        assert compare_files(
            OUT_FILE_DEFAULT, expected_file
        ), "Output file does not match expected output file"

    """
    Run the local neighborhood algorithm with a missing input file
    """

    def test_missing_file(self):
        with pytest.raises(OSError):
            LocalNeighborhood.run(
                nodes=TEST_DIR+'input/ln-nodes.txt',
                network=TEST_DIR+'input/missing.txt',
                output_file=OUT_FILE_DEFAULT
            )

    """
    Run the local neighborhood algorithm with an improperly formatted network file
    """

    def test_format_error(self):
        with pytest.raises(ValueError):
            LocalNeighborhood.run(
                network=TEST_DIR+'input/ln-bad-network.txt',
                output_file=OUT_FILE_DEFAULT
            )

    # Write tests for the Local Neighborhood run function here

    @pytest.mark.skipif(not shutil.which('singularity'), reason='Singularity not found on system')
    def test_ln_required(self):
        out_path = Path(OUT_FILE_DEFAULT)
        out_path.unlink(missing_ok=True)
        # Only include required arguments
        LocalNeighborhood.run(
            nodes=TEST_DIR+'input/ln-nodes.txt',
            network=TEST_DIR+'input/ln-network.txt',
            output_file=OUT_FILE_DEFAULT
        )
        assert out_path.exists()
