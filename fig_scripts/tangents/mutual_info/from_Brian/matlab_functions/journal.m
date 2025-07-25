function Need2Read = journal(Num)
%Need2Read = journal(Num)

NumDays = datenum(date)-datenum(2006,7,20);
Need2Read = (NumDays*5/7 - Num);