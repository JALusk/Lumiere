test:
	python -Wall -m unittest discover -v

clean:
	rm -rf build/
	rm -rf *.egg-info
	rm -rf dist/
