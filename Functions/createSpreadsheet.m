%%
% NAME: CREATE SPREADSHEET
% AUTHOR: JANE HUMPHREY (janehumphrey@outlook.com)

function createSpreadsheet(excel,file,text,footerRows)

if nargin<3
    error('Not enough input arguments.');
elseif nargin==3
    footerRows = 0;
end

[nRows,nColumns] = size(text);
nRowsStr = num2str(nRows);
if nColumns<=26
    nColumnsStr = char('A'-1+nColumns);
elseif nColumns<=26*27
    part1 = char('A'-1+floor((nColumns-1)/26));
    part2 = char('A'+mod(nColumns-1,26));
    nColumnsStr = [part1,part2];
else
    part1 = char('A'-1+floor((nColumns-1)/(26*27)));
    part2 = char('A'+floor(mod(nColumns-1,26*27)/26));
    part3 = char('A'+mod(nColumns-1,26));
    nColumnsStr = [part1,part2,part3];
end

workbook = excel.Workbook.Add;
sheet = workbook.ActiveSheet;
sheet.Activate;

rangeStr = ['A1:',nColumnsStr,nRowsStr];
range = sheet.Range(rangeStr);
range.Value = text;
headerStr = ['A1:',nColumnsStr,'1'];
header = sheet.Range(headerStr);
header.Font.Bold = 1;
if footerRows>1
    firstRow = nRows-(footerRows-1);
    firstRowStr = num2str(firstRow);
    footerStr = ['A',firstRowStr,':',nColumnsStr,nRowsStr];
    footer = sheet.Range(footerStr);
    footer.Font.Bold = 1;
end

sheet.Rows.VerticalAlignment = -4108;
sheet.Rows.AutoFit;
sheet.Columns.AutoFit;
workbook.SaveAs(file);
workbook.Close;