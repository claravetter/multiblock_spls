%% labels hammers

% 
% for i=1:size(labelsdartelhammers,1)
%     temp = labelsdartelhammers{i};
%     labels_hammers{i,1} = temp((strfind(temp, '<label><index>')+size('<label><index>',2)):(strfind(temp, '</index>')-1));
%     labels_hammers{i,2} = temp((strfind(temp, 'index><short_name>')+size('index><short_name>',2)):(strfind(temp, '</short_name><name>')-1));
%     labels_hammers{i,3} = temp((strfind(temp, '/short_name><name>')+size('/short_name><name>',2)):(strfind(temp, '</name><RGBA')-1));
%     labels_hammers{i,4} = temp((strfind(temp, 'XYZmm>')+size('XYZmm>',2)):(strfind(temp, '</XYZmm></label>')-1));
% end

