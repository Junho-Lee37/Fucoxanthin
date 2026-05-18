function genes = read_gene_list(fname, G)
%READ_GENE_LIST read genes_used.txt (1 gene per line) in same order as params
genes = cell(G,1);

if exist(fname,'file')
    txt = readlines(fname);
    txt = txt(txt ~= "");
    n = min(numel(txt), G);
    for i=1:n
        genes{i} = char(strtrim(txt(i)));
    end
    for i=n+1:G
        genes{i} = sprintf('gene_%d', i);
    end
else
    warning('genes file not found: %s. Using generic names.', fname);
    for i=1:G
        genes{i} = sprintf('gene_%d', i);
    end
end
end
