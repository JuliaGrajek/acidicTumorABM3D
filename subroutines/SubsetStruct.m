function T = SubsetStruct(S, idx)
    fields = fieldnames(S);
    for i = 1:numel(fields)
        Field    = fields{i};
        T.(Field) = S.(Field)(idx);
    end
end