%% function to check for inverse weightings in opt parameters vectors

function [output_matrix, inv] = dp_inverse_array(input_matrix, inv, parameters_names)

output_matrix=input_matrix;

if ~isfield(inv, 'log')
    
    if ~isfield(inv, 'to_extract')
        disp('Please choose which feature to check: u or v?');
    else
        n = strcmp(parameters_names,inv.to_extract);
    end
    
    temp_extract = input_matrix(:,n);

    val=[]; ind=[];
    for i=1:size(temp_extract,1)
        [val(i,1) ind(i,1)] = max(abs(temp_extract{i}));
    end
    
    if range(ind)==0
        decision_numbers = ind(1,1);
    else
        decision_numbers = (unique(ind));
    end
    
    
    temp=[];
    for i=1:size(input_matrix,1)
        temp{i,1} = (temp_extract{i,:}(decision_numbers) > 0)*1;
        temp{i,2} = (temp_extract{i,:}(decision_numbers) < 0)*(-1);
    end
    
    pos = sum([temp{:,1}],2) > 5;
    neg = sum([temp{:,2}],2) < -5;
    
    inv.log=true;
    for i=1:size(temp_extract,1)
        decision_1 = sum(temp_extract{i,:}(decision_numbers(pos))<0) > 0;
        decision_2 = sum(temp_extract{i,:}(decision_numbers(neg))>0) > 0;
        if decision_1==1 || decision_2==1
            output_matrix{i,n} = (-1)*temp_extract{i};
            inv.log(i,1) = true;
            %     elseif decision_1==0 && decision_2 ==0
            %         output_matrix{i,n} = (-1)*temp_extract{i};
            %         inv_log(i,1) = true;
        else
            output_matrix{i,n} = temp_extract{i};
            inv.log(i,1) = false;
            %         disp(['There is a problem with opt parameter ' num2str(i) '. Please check it manually']);
        end
    end
    
else
    
    if ~isfield(inv, 'to_apply')
        disp('Please choose which feature to check: u or v?');
    else
        n = strcmp(parameters_names,inv.to_apply);
    end
    
    temp_extract = input_matrix(:,n);
    
    for i=1:size(inv.log,1)
        if inv.log(i)
            output_matrix{i,n} = (-1)*temp_extract{i};
        end
    end
end

end