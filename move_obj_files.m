function dwi_object=move_obj_files(dwi_object,folder_to_write,name)
%right now I'm gonna implement it as if every element has the same fn. we
%can complicate later

% really this is happening bc in most cases, each dwi obj element doesn't
% need it's own fn field (like they dont need their own json fields), but
% sometimes we concatenate dwi objs from different sources and we need to
% handle these cases. The thing that makes this all a little silly is that
% there is minimal handling of these cases -- yet.

fn_field_names=fieldnames(dwi_object(1).fns);
for i=1:numel(fn_field_names) % for each type of file
    allfns=cell(numel(dwi_object),1); %initialize cell array
    for j=1:numel(allfns)
        if isfield(dwi_object(j).fns,fn_field_names{i})
            allfns{j}=dwi_object(j).fns.(fn_field_names{i}); % stick the file name in this array
        end
    end
    
    %handling weird input conditions -- might need some changing
    % best idea is to move each unique fns, make sure only to try on
    % nonempty fn names. but where do you move it to?
    if ~isequal(allfns{:})
        warning('your filenames are not all the same for %s, using only the first',fn_field_names{i})
    end
    allfnsempty=cellfun(@(x) isempty(x),allfns);
    allfns=allfns(~allfnsempty);    
    
    %move
    fprintf('moving %s\n',fn_field_names{i})
    mkdir(folder_to_write)
    new_fn=fullfile(folder_to_write,[name '.' fn_field_names{i}]);
    movefile(allfns{1},new_fn)
    
    %change dwi_object to reflect new location
    dwi_object=arrayfun(@(x) setfield(x,'fns',fn_field_names{i},new_fn),dwi_object);

end
    