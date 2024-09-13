function grid = create_grid_define_range(density, varargin)
    % Set default values for start_value and stop_value
    start_value = 1;
    stop_value = 0;
    
    % Parse varargin for named arguments
    i = 1;
    while i <= length(varargin)
        switch varargin{i}
            case 'start_value'
                if i + 1 <= length(varargin)
                    start_value = varargin{i + 1};
                    i = i + 2;
                else
                    error('Missing value for start_value.');
                end
            case 'stop_value'
                if i + 1 <= length(varargin)
                    stop_value = varargin{i + 1};
                    i = i + 2;
                else
                    error('Missing value for stop_value.');
                end
            otherwise
                error('Unknown parameter name: %s', varargin{i});
        end
    end
    
    % Create the grid struct
    grid = struct('start', start_value, 'end', stop_value, 'density', density); 
end
