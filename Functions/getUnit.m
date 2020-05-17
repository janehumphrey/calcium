%%
% NAME: GET UNIT
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function [unitStr] = getUnit(unit)

microStr = char(hex2dec('03BC'));
minusStr = char(hex2dec('207B'));
oneStr = char(hex2dec('00B9'));
twoStr = char(hex2dec('00B2'));
spaceStr = char(hex2dec('2009'));

if strcmp(unit,'ug')==1
    unitStr = [microStr,'g'];
elseif strcmp(unit,'um')==1
    unitStr = [microStr,'m'];
elseif strcmp(unit,'um2')==1
    unitStr = [microStr,'m',twoStr];
elseif strcmp(unit,'perS')==1
    unitStr = ['s',minusStr,oneStr];
elseif strcmp(unit,'umPerS')==1
    unitStr = [microStr,'m',spaceStr,'s',minusStr,oneStr];
elseif strcmp(unit,'um2PerS')
    unitStr = [microStr,'m',twoStr,spaceStr,'s',minusStr,oneStr];
elseif strcmp(unit,'molPerUm2')
    unitStr = ['molecules',spaceStr,microStr,'m',minusStr,twoStr];
else
    error('Unknown unit.');
end