function parFormatDescription
% Description of .par file format for specifying event-related or block-
% design stimuli.
%
% The format is as follows: a tab-delimited ASCII text file with the
% following columns:
%
% (1) [\tab] (2) [\tab] (3) [\tab] (4) [\tab] (5) [\tab] (6) [\tab] (7)
%
% where: (1) is the onset time in seconds;
%        (2) is an integer specifying an event type, with 0 specifying
%            a 'baseline' or 'null' condition;
%        (3) is a text label specifying the name of a condition. (Only
%            the first occurrence of a condition is parsed to get the 
%            name);
%        (4) specifies a color to associate w/ each condition, for use
%            in visualizing data that is sorted by condition. The color
%            is specified in matlab format, either a single letter
%            (e.g. 'w' for 'white', 'r' for 'red', 'm' for 'magenta')
%            or the format '[R G B]' where R, G, and B are the relative
%            intensities of the red, green and blue channels, each ranging
%            from 0-1 (see 'help plot');
%        (5) specifies a path to an image to associate w/ each condition,
%            for use in visualizing data that is sorted by condition;
%        (6) and (7) specify behavioral data, such as response times
%            or correct/incorrect responses.
%
%   All columns except for (1) and (2) are optional. If you don't have
%   information for a column, but would like to include it for a later
%   column, you can add the text '' or [] for that column (or add two
%   tabs).
%
%
%
% Importantly, .par files here have an implicit assumption that 
% there's always _something_ going on during an experiment: events don't 
% have a separately-specified duration, but are assumed to continue
% until the next event starts. Therefore, you may want to specify
% explicitly the end of a scan as a separate event (which can be 
% assigned to a null condition or a negative number, if you like).

% if someone called this as a function, just give them the info:
p = fileparts(mfilename);
web(fullfile(p,'parFormatDescription.txt'));

return