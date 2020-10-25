VERSION=$(shell python3 -c "from configparser import ConfigParser; p = ConfigParser(); p.read('setup.cfg'); print(p['metadata']['version'])")

default:
	@echo "\"make publish\"?"

tag:

	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi
	curl -H "Authorization: token `cat $(HOME)/.github-access-token`" -d '{"tag_name": "v$(VERSION)"}' https://api.github.com/repos/krober10nd/SeismicMesh/releases

upload:
	# Make sure we're on the master branch
	@if [ "$(shell git rev-parse --abbrev-ref HEAD)" != "master" ]; then exit 1; fi

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
	isort -rc SeimicMesh/ tests/*.py
	black SeismicMesh/ benchmarks/*.py tests/*.py
	clang-format -i SeismicMesh/generation/cpp/*.cpp SeismicMesh/migration/cpp/*.cpp SeismicMesh/sizing/cpp/*.cpp SeismicMesh/geometry/cpp/*.cpp

black:
	black .

lint:
	isort --check . SeismicMesh/ setup.py benchmarks/*.py tests/*.py
	black --check SeismicMesh/ setup.py benchmarks/*.py tests/*.py
	flake8 setup.py SeismicMesh/ benchmarks/*.py tests/*.py
