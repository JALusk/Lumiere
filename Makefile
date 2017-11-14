test:
	python3.6 -m unittest discover -s tests.unit -v

integrate:
	python3.6 -m unittest discover -s tests.integration -v
