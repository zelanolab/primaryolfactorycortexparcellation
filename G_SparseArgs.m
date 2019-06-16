function options = G_SparseArgs( options, varargin)
% function sparse key-value pairs input
% 
% This function only validates user input Name-Value pairs
% 
% Input
%   options: structure with keys and default values
% 	args: user input key-value pairs
% 
% Output
% 	options: a structure of which each field has a user-specified value
% 
% Example Usage
%   	options = struct( 'seg_len', [],...
%                          'seg_overlap', []);
%       args = {'seg_len', 1024, 'seg_overlap', 512}; % Matlab varargin
%       options = G_SparseArgs( options,  args)
%
% 

%% validate input arguments
if nargin < 2
    error( 'Not enough input argument. Use as G_SparseArgs( options,  args);');    
elseif nargin > 2
    error( 'Too many input arguments. Use as G_SparseArgs( options,  args);');    
else
    % do nth
end

if ~isstruct( options)
    error( 'options must be a structure.');
end

%% sparse input arguments
option_names = fieldnames( options);

args = varargin{1};
nbargus = length( args);
if mod( nbargus, 2) ~= 0
   error('Key-Value pairs needed.')
end

for pair = reshape( args, 2, []) 
   in_name = lower( pair{1}); 
   if any( strcmpi( in_name, option_names))
      options.( in_name) = pair{2};
   else
      error('%s is not a recognized parameter name.', in_name)
   end   
end

end % function