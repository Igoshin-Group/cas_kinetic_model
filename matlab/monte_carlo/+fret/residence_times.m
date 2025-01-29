function T = residence_times(trajectory, num_states)
T = cell(num_states,1);

t = trajectory.time;
S = trajectory.state;

if isempty(t)
    return;
end

% find state residence times
prev_state = S(1);
entry_time = t(1);
for ii = 2:length(t)
    cur_state = S(ii);
    if (cur_state ~= prev_state)
        % compute the total residence time
        exit_time = t(ii);
        residence_time = exit_time - entry_time;
        
        T{prev_state}(end+1) = residence_time;
        
        % update the entry time and previous state
        prev_state = cur_state;
        entry_time = exit_time;
    end
end
end