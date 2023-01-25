function [source_train, target_train, source_test, target_test] = ...
    embed4D(source_data, target_data, train_span, test_span, embed_step, embed_before, time_step, nahead)
% this version is for 4D data.
%
% includes nahead as input.

embed_size = embed_before + 1;
M = size(source_data, 2);  % number of time series used for source signal

% set embedded data: source_train, target_train
% columns for source_train: data1(t-embed_before*embed_step),,,,,, data1(t),data2(t-embed_before*embed_step),,,,,, data2(t),

start_ind = train_span(1) + embed_before * embed_step;
end_ind = train_span(2) - nahead;
%disp([start_ind, end_ind])

ind_seq = start_ind:time_step:end_ind; % sequence of time indices to pick up for training
source_train = NaN(length(ind_seq), embed_size*M);
columns = 1:embed_size:(1+embed_size*(M-1));  % columns of source_train to put each timeseries
% embed
for shift = ((-embed_before):0)*embed_step
    %disp(['shift = ', num2str(shift)])
    source_train(:, columns) = source_data(ind_seq + shift, :);   % from old to now
    columns = columns + 1; % add to each column number
end
target_train = target_data(ind_seq+nahead,:);


% set embedded data, source_test, target_test
% columns for target_test: same as above
if length(test_span) < 2
    source_test = [];
    target_test = [];
else

    start_ind = test_span(1) + embed_before * embed_step;
    end_ind = test_span(2) - nahead;
    %disp([start_ind, end_ind])

    ind_seq = start_ind:time_step:end_ind; % sequence of time indices to pick up
    source_test = NaN(length(ind_seq), embed_size*M);
    columns = 1:embed_size:(1+embed_size*(M-1));
    % embed
    for shift = ((-embed_before):0)*embed_step
        %disp(['shift = ', num2str(shift)])
        source_test(:, columns) = source_data(ind_seq + shift, :);
        columns = columns + 1;
    end
    target_test = target_data(ind_seq+nahead,:);

end

end
