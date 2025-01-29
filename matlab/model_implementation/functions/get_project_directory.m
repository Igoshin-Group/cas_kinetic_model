function dir = get_project_directory()
cd(".."); % jump to code directory
dir = pwd(); % get the working directory
end