function simplewarning(msg)
% Print orange text without adding stack trace and line numbers
% Based upon the following hack:
% https://undocumentedmatlab.com/articles/another-command-window-text-color-hack

msg = sprintf('Warning: %s\n]',msg);
fprintf(['[', 8, msg, 8]);
end