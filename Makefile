test:
	python3 -m unittest discover -s tests.unit -v

integrate:
	python3 -m unittest discover -s tests.integration -v

parallelwiggle:
	mpiexec python3 -m unittest tests.unit.test_calc_wiggled -v

