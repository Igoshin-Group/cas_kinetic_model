function distribution = determine_transition_path_distribution(trajectory)
    time = trajectory.time;
    state_trajectory = trajectory.state;
    trajectory_length = length(state_trajectory);
    fret_state = state_trajectory;
    distribution = [];
    
    entry_time = 0;
    exit_time = 0;
    has_exited = false;

    for ii = 2:trajectory_length
        state = fret_state(ii);
        prev_state = fret_state(ii-1);
        
        if (state == 1)
            has_exited = false;
            continue;
        end
        
        if ((prev_state == 1) && (state == 2))
            entry_time = time(ii);
            has_exited = false;
            continue;
        end
        
        if (state == 2)
            continue;
        end
        
        if ((prev_state == 2) && (state == 3) && (~has_exited))
            exit_time = time(ii); % need to have time to state 3
            distribution(end+1) = exit_time - entry_time;
            has_exited = true;
            continue;
        end
    end
end