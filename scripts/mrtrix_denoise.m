function varargout=mrtrix_denoise(fn_4d,suffix,shell,exe)
% operates on files/folders
%
% This function is included for ease of use on windows -- dwidenoise is much
% more straightforward to use on unix-like systems. Just use command line for these cases

%convert strings to char just in case
if isstring(fn_4d)
    fn_4d=char(fn_4d);
end
if isstring(suffix)
    suffix=char(suffix);
end

%this only works on 4d volumes, so:
assert(numel(spm_vol(fn_4d))>1,'only works on 4d niis');

%% do the simple case if shell,exe not defined
%note: I wrote this in vim and not my normal matlab environment, so be wary of things like unclosed quotes etc
if ~exist('shell','var') && ~exist('exe','var') 
    warning('no shell or exe given -- assuming you''re on a mac or unix and you''ve got mrtrix set up for command line')
    [a1,b1,c1]=fileparts(fn_4d);

    %commands
    command1=['dwidenoise "' fullfile(a1,[b1 c1]) '" "' fullfile(a1,[b1 suffix c1]) '" -noise "' fullfile(a1,['noise.'...
    'nii']) '" -force'];
    [status{1},cmdout{1}] = system(command1);

    command2= ['mrcalc "' fullfile(a1,[b1 c1]) '" "' fullfile(a1,[b1 suffix c1]) '" -subtract'...
    ' "' fullfile(a1,'res.nii') '" -force'];
    [status{2},cmdout{2}] = system(command2);

    %handle outs -- wrote this a while ago, might not be necessary now
    if nargout==2
        varargout{1}=status;
        varargout{2}=cmdout;
    elseif nargout==1
        varargout{1}=cmdout;
    end
    %that's all folks!
    return
end
%convert strings to char just in case

if isstring(shell)
    shell=char(shell);
end
if isstring(exe)
    exe=char(exe);
end

%find out if the user has given you the path to dwidenoise or to the bin
%folder
dd=dir(exe);
if dd(1).isdir %if it IS a directory
    casee=1;
elseif ~dd(1).isdir %it's not a directory; assuming it's the actual exe itself
    casee=2;
elseif isempty(dd(1))
    error('dir call to exe came up empty')
else
    error('exe path isn''t parsed correctly, talk to andrew')
end

[a1,b1,c1]=fileparts(fn_4d);
[a2,b2,c2]=fileparts(exe); % grab mrcalc with the same path as dwidenoise; this should work in almost every case

%command1 = ['C:/Programs/msys64/mingw64.exe C:\Programs\msys64\home\Andrew\mrtrix3\bin\dwidenoise.exe ' root '/com/denoising/4dall.nii ' root '/com/denoising/4D_DN.nii -noise ' root '/com/denoising/noise.nii'];
%[status,cmdout] = system(command1);

if isempty(shell)
    sep='';
else
    sep=' ';
end
if casee==2 %exe dwidenoise is given
    command1=[shell sep exe ' "' fullfile(a1,[b1 c1]) '" "' fullfile(a1,[b1 suffix c1]) '" -noise "' fullfile(a1,'noise.nii') '" -force']
    [status{1},cmdout{1}] = system(command1);
end

if casee==1 %exe folder is given
    command1=[shell sep fullfile(exe,'dwidenoise') ' "' fullfile(a1,[b1 c1]) '" "' fullfile(a1,[b1 suffix c1]) '" -noise "' fullfile(a1,'noise.nii') '" -force']
    [status{1},cmdout{1}] = system(command1);
end
%make sure fullfile doesn't screw this up; the file/directory seperation
%character that fullfile uses could be wrong, because msys shell uses
%different character than OS.
disp('MATLAB continues after calling EXE')
% Check to see if the EXE process exists
flag = true;
while flag
     disp('EXE still running');   
     flag = isprocess('mintty.exe');
     pause(1) 
end
disp('EXE Done')
disp('continuing with MATLAB script')

fnarr
%command = ['C:/Programs/msys64/mingw64.exe C:\Programs\msys64\home\Andrew\mrtrix3\bin\mrcalc ' root '/com/denoising/4dall.nii ' root '/com/denoising/4D_DN.nii -subtract ' root '/com/denoising/res.nii -force'];
%[status,cmdout] = system(command);
if casee==2 %exe is given directly
    command2= [shell sep fullfile(a2,'mrcalc') ' "' fullfile(a1,[b1 c1]) '" "' fullfile(a1,[b1 suffix c1]) '" -subtract "' fullfile(a1,'res.nii') '" -force'];
    [status{2},cmdout{2}] = system(command2);
end
if casee==1 %exe folder is given
    command2= [shell sep fullfile(a1,'mrcalc') ' "' fullfile(a1,[b1 c1]) '" "' fullfile(a1,[b1 suffix c1]) '" -subtract "' fullfile(a1,'res.nii') '" -force'];
    [status{2},cmdout{2}] = system(command2);
end

disp('MATLAB continues after calling EXE')
% Check to see if the EXE process exists
flag = true;
while flag
     disp('EXE still running');   
     flag = isprocess('mintty.exe');
     pause(1) 
end
disp('EXE Done')
disp('continuing with MATLAB script')


if nargout==2
    varargout{1}=status;
    varargout{2}=cmdout;
elseif nargout==1
    varargout{1}=cmdout;
end
end
