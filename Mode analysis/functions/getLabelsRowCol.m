function [Label_row, Label_col] = getLabelsRowCol(Q_temp)
    thresh = 30; % nm

    if length(size(Q_temp))==3
        Q_temp = mean(Q_temp,3);
    end
    Label_row = zeros(size(Q_temp,1),1);
    Label_col = zeros(size(Q_temp,1),1);
    % row
    [a,b] = sort(Q_temp(:,2));
    for i = 1:length(a)
        if i == 1
            Label_row(b(i)) = 1;
            continue
        end
        if a(i)-a(i-1) > thresh
            Label_row(b(i)) = Label_row(b(i-1))+1;
        else
            Label_row(b(i)) = Label_row(b(i-1));
        end
    end
    % col
    rows = unique(Label_row);
    for i = rows'
        select = find(Label_row == i);
        temp = Q_temp(select,1);
        [~,b] = sort(temp);
        for j = 1:length(b)
            Label_col(select(b(j))) = j;
        end
    end

    % figure(); hold on;
    % scatter(Q_temp(:,1),Q_temp(:,2))
    % for i = 1:max(Label_col)
    %     for j = 1:max(Label_row)
    %         select = find(Label_row == j & Label_col == i);
    %         text(Q_temp(select,1),Q_temp(select,2),[num2str(i),',',num2str(j)])
    %     end
    % end
end