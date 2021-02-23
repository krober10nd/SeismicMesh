import pytest

from SeismicMesh import (
    Rectangle,
    Disk,
    Ball,
    Cube,
    Torus,
    Prism,
    Cylinder,
    generate_mesh,
)


def test_rotations2d():
    bbox = (0.0, 1.0, 0.0, 1.0)
    geo = []
    geo.append(Disk((0.5, 0.5), 0.1, rotate=0.2 * 3.14))
    geo.append(Rectangle(bbox, rotate=0.3 * 3.14))

    for domain in geo:
        p, c = generate_mesh(domain=domain, edge_length=0.1)


def test_rotations3d():
    bbox = (0.0, 1.0, 0.0, 1.0, 0.0, 1.0)
    geo = []
    # geo.append(Ball((0.5, 0.5, 0.5), 0.1, rotate=0.0 * 3.14))
    geo.append(Cube(bbox, rotate=0.4 * 3.14))
    # geo.append(Torus(1.0, 0.5, rotate=0.4 * 3.14))
    # geo.append(Prism(0.5, 0.5, rotate=0.5 * 3.14))
    # geo.append(Cylinder(1.0, 0.5, rotate=0.6 * 3.14))

    for domain in geo:
        p, c = generate_mesh(domain=domain, edge_length=0.1)


@pytest.mark.serial
def test_rotations():
    test_rotations2d()
    test_rotations3d()


if __name__ == "__main__":
    test_rotations()
