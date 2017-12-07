function charindvars = symbolArrayToString(indvars)

indvarsAsStrs = arrayfun(@(arg) char(arg), indvars, 'UniformOutput', false);
charindvars = strcat('[', strjoin(indvarsAsStrs, ','), ']');

end

