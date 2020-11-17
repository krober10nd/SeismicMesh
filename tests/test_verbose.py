import os
import sys

import pytest
import SeismicMesh


@pytest.mark.serial
def test_verbose():
    square = SeismicMesh.Rectangle((0.0, 1.0, 0.0, 1.0))

    for verbosity, correct_size in zip([0, 1, 2], [0, 246, 6068]):
        sys.stdout = open("output.txt", "w")
        points, cells = SeismicMesh.generate_mesh(
            domain=square, edge_length=0.1, verbose=verbosity
        )
        sys.stdout.close()
        sys.stdout = sys.__stdout__
        output_size = os.path.getsize("output.txt")

        assert output_size == correct_size


if __name__ == "__main__":
    test_verbose()
