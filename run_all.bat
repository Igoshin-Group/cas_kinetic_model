@echo off
echo Running FRET analysis.
call run_fret_analysis.bat
echo Running Monte Carlo simulations.
call run_monte_carlo_simulations.bat
echo Running model analysis.
call run_model_analysis.bat
