test:
	python3 -m unittest discover -s tests.unit -v

integrate:
	python3 -m unittest discover -s tests.integration -v
