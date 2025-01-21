% Put this code in startup.m in your Matlab home folder to execute it on
% each Matlab launch. This folder can be found by typing 'userpath' in your
% console. If this file does not exist you can manually create it.
%
% Another way to use it is to put 'better_colorwheel' in your plotting
% scripts.
%
% Sebastiaan Goossens
% November 2019

colors = [

% Primary color set
031 119 180
214 039 040
255 127 014
044 160 044
148 103 198
140 086 075
227 119 194
023 190 207
127 127 127
188 189 034

% Secondary color set
% Lighter version of primary
174 199 232
255 152 150
255 187 120
152 223 138
197 176 213
196 156 148
247 182 210
158 218 229
199 199 199
219 219 141

% Matlab wants floating point numbers between 0 and 1
]/255;

% Set default colors
set(groot,'defaultAxesColorOrder',colors);