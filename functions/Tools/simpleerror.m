function simpleerror(msg)
% Print red text without adding stack trace and line numbers
% Based upon the following hack:
% https://nl.mathworks.com/matlabcentral/answers/15328-adding-color-in-command-window-output

msg = sprintf('Error: %s\n',msg);
fprintf(2, msg);
end