% Supports heat map plots for the file Experiment4.m
% See Experiment4.m for more details
% Edits the correlation structure matrices of the format in [1]
function data_ed = data_edit_index(data)
data_ed = zeros(size(data));
for i=1:size(data,1)
    j=size(data,1) - i+1;
    data_ed(j,:) = data(i,:);
end
end