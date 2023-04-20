test:
	python3 -m unittest discover -s tests.unit -v

integrate:
	python3 -m unittest discover -s tests.integration -v

parallelwiggle:
	mpiexec -n 5 python3 -m superbol.mpi_shell -v

