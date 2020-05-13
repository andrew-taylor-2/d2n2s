function copy_obj_files(dwi_object,folder_to_write,name)
%right now I'm gonna implement it as if every element has the same fn. we
%can complicate later
fn_field_names=fieldnames(dwi_object(1).fns);
for i=1:numel(fn_field_names)
    allfns=cell(numel(dwi_object),1);
    for j=1:numel(allfns)
        if isfield(dwi_object(j).fns,fn_field_names{i})
            allfns{j}=dwi_object(j).fns.(fn_field_names{i});
        end
    end
    if ~isequal(allfns{:})
        warning('your filenames are not all the same for %s, using only the first',fn_field_names{i})
    end
    allfnsempty=cellfun(@(x) isempty(x),allfns);
    allfns=allfns(~allfnsempty);    
    
    fprintf('copying %s\n',fn_field_names{i})
    mkdir(folder_to_write)
    copyfile(allfns{1},fullfile(folder_to_write,[name '.' fn_field_names{i}]))
end
    