VERSION=$(shell python3 -c "from configparser import ConfigParser; p = ConfigParser(); p.read('setup.cfg'); print(p['metadata']['version'])")

default:
	@echo "\"make publish\"?"

tag:

	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "par3d" ]; then exit 1; fi
	curl -H "Authorization: token `cat $(HOME)/.github-access-token`" -d '{"tag_name": "v$(VERSION)"}' https://api.github.com/repos/krober10nd/SeismicMesh/releases

upload:
	# Make sure we're on the par3d branch
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "par3d" ]; then exit 1; fi

	rm -rf dist/*
	python3 setup.py sdist bdist_wheel
	twine upload dist/*

publish: tag upload

clean:
	@find . | grep -E "(__pycache__|\.pyc|\.pyo$\)" | xargs rm -rf
	@rm -rf build/*
	@rm -rf SeismicMesh.egg-info/
	@rm -rf dist/

format:
	isort -rc SeimicMesh/ examples/*.py tests/*.py
	black SeismicMesh/ examples/*.py tests/*.py
	blacken-docs README.md
	clang-format -i SeismicMesh/generation/cpp/*.cpp SeismicMesh/migration/cpp/*.cpp SeismicMesh/sizing/cpp/*.cpp SeismicMesh/geometry/cpp/*.cpp

black:
	black .

lint:
	isort --check . SeismicMesh/ setup.py examples/*.py tests/*.py
	black --check SeismicMesh/ setup.py examples/*.py tests/*.py
	flake8 setup.py SeismicMesh/ examples/*.py tests/*.py
