@ECHO OFF

:test
python -m unittest discover -s tests.unit -v
goto end

:integrate
python -m unittest discover -s tests.integration -v
goto end

:end