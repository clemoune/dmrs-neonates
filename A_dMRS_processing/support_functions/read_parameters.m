function parameters = read_parameters(filename)

%
% This function is part of:
%
% Bruker - Graphical user interface to image Bruker data 
%
% Usage - call bruker from Matlab console
%
% Written by Matteo Caffini, PhD
% Dipartimento di Elettronica, Informatica e Bioingegneria
% Politecnico di Milano, Milano, ITALY
%
% Copyright (C) 2014 Matteo Caffini <matteo.caffini@polimi.it>
%

parameters = [];

% Prompt for a file if not given as an input argument
if nargin == 0
    [fn,fp] = uigetfile({'*.*','All Files (*.*)'},'Open a Bruker parameters file (method, acqp, ...)');
    if isequal(fn,0)
        return
    end
    filename = [fp,fn];
elseif nargin > 1
    error('Too many input arguments.');
end

% Open the file for reading
fileId = fopen(filename,'r');
if fileId < 0
    error('Could not open file "%s" for reading.',filename);
end

% Check that the file is a Bruker parameters file
str = fread(fileId,20,'char=>char');
if isempty(regexp(str.','^\s*##TITLE'))
    fclose(fileId);
    error('File "%s" is not a valid Bruker parameters file.',filename)
end
fseek(fileId,0,-1); % Rewind file

C = fread(fileId,inf,'char');
fclose(fileId);

% Remove carriage returns
C(C==13)=[];

% Convert to string
C = char(C.');

% Remove comment lines
C = regexprep(C,'\$\$([^\n]*)\n','');

% Remove unnecessary line breaks
f = @remove_line_breaks;
C=regexprep(C,'^(\s*[^#].*?)(?=\n\s*#)','${f($1)}','lineanchors');
C=regexprep(C,'(\([^\)]+?)\n(.*?\))','${f([$1,$2])}','lineanchors');
CC = regexp(C,'\s*##','split');
CC(1)=[];

% Parse the file line-by-line
for ii=1:length(CC)
    
    str = CC{ii};
    if strncmp(str,'END=',4)
        continue
    end
    
    % The commented regexp sometimes fails with long strings...
    %param = regexp(str,'^(.*)=','tokens','once');
    ind = find(str==61); % Find '=' chars...
    if isempty(ind)
        param='';
    else
        param=str(1:ind(1)-1);
    end
    %param = strrep(param{1},'$','');
    param = strrep(param,'$','');
    param = check_parameter_string(param);
    
    if any(str==sprintf('\n'))
        % Get size
        sz = regexp(str,'=\s*\((.*?)\)\s*\n','tokens','once');
        sz = str2num(['[',sz{1},']']);
        
        % Parse value
        value = regexp(str,'\n(.*)$','tokens','once');
        value = value{1};
        value = check_parameter_value(value,sz);
    else
        value = regexp(str,'=\s*(.*)','tokens','once');
        value = value{1};
        value = check_parameter_value(value);
    end
    
    % Add to structure
    parameters.(param) = value;
    
end

%% Nested functions

    function out = remove_line_breaks(str)
        out = strrep(str,sprintf('\n'),'');
    end

    function out = check_parameter_value(val,sz)
        
        if nargin == 1
            sz = 0;
        end
        
        % Remove insignificant whitespace
        val = strtrim(val);
        
        if isempty(val)
            out = val;
            return
        end
        
        % Handle strings and string lists
        if val(1) == '<' && val(end) == '>'
            val(val=='<')='''';
            val(val=='>')='''';
            out = eval(['{',val,'}']);
            if length(out) == 1
                out = out{1};
            end
            return
        end
        
        % Handle cell matrices
        if val(1) == '(' && val(end) == ')'
            nRows = length(find(val==')'));
            
            % Nested tables are not supported. This is a workaround for nested tables
            % and everything is read in a single lined table...
            if nRows ~= sz && sz>0
                nRows=sz;
            end
            
            val(1) = '';
            val(end) = '';
            val(val=='(')='';
            val(val==')')=',';
            val(val=='<')='';
            val(val=='>')='';
            
            % Split using the commas
            val_split = regexp(val,',\s+','split');
            val_out = cell(size(val_split));
            
            % Try to convert to numbers
            for ii = 1:length(val_split)
                num = str2double(val_split{ii});
                if isnan(num)
                    val_out{ii} = val_split{ii};
                else
                    val_out{ii} = num;
                end
            end
            
            
            out = reshape(val_out,[],nRows).';
            return
        end
        
        % Check if the string contains only numbers before tryin to convert to a
        % number. str2num uses eval command and if the string matches to a
        % function name strange things can happen...
        tmp2 = regexp(val,'[^\d\.\seE-+]');
        if ~isempty(tmp2)
            out = val;
            return
        end
        
        % Convert value to numeric if possible
        tmp = str2num(val);
        if ~isempty(tmp) && isreal(tmp)
            if length(sz)>1
                tmp = reshape(tmp,sz(2),sz(1),[]);
                tmp = permute(tmp,[2 1 3]);
            end
            out = tmp;
            return
        end
        
        out = val;
    end

    function out = check_parameter_string(param)
        
        alphabets = 'abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ';
        numbers = '1234567890';
        
        % Remove insignificant whitespace
        param = strtrim(param);
        
        if isempty(param)
            out = 'EMPTY_PARAM';
            return
        end
        
        % Check parameter starts with a valid structure field character
        if ~any(param(1)==alphabets)
            param = ['PAR_',param];
        end
        
        % Check that the parameter string does not contain any illegal characters
        % (for Matlab structure fields)
        ind = ~ismember(param,[alphabets,numbers,'_']);
        if any(ind)
            param(ind) = '_';
        end
        
        out = param;
    end

end