test:
	python3.6 -Wall -m unittest discover -s tests.unit -v

integrate:
	python3.6 -Wall -m unittest discover -s tests.integration -v
