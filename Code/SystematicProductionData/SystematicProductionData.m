% Data obtained from: http://www.spacelaunchreport.com/logyear.html
totalLaunchData = [
        2000  78  7  85;...
        2001  51  8  59;...
        2002  58  7  65;...
        2003  59  4  63;...
        2004  52  2  54;...
        2005  51  4  55;...
        2006  61  5  66;...
        2007  63  5  68;...
        2008  61  7  68;...
        2009  69  9  78;...
        2010  67  7  74;...
        2011  77  7  84;...
        2012  73  5  78;...
        2013  76  5  81;...
        2014  92  4  92;...
        2015  86  4  86];
    
LEOLaunchData = [
        2000  33;...
        2001  27;...
        2002  31;...
        2003  25;...
        2004  26;...
        2005  29;...
        2006  33;...
        2007  36;...
        2008  36;...
        2009  45;...
        2010  37;...
        2011  43;...
        2012  40;...
        2013  48;...
        2014  50;...
        2015  44];  
    
LEOISSData=[
        2010  12;...
        2011  13;...
        2012  12;...
        2013  12;...
        2014  13;...
        2015  12];  
LEOSData=[
        2010  12;...
        2011  17;...
        2012  14;...
        2013  16;...
        2014  23;...
        2015  18];  
    
LEOOtherData=[
        2010  14;...
        2011  13;...
        2012  14;...
        2013  20;...
        2014  14;...
        2015  11]; 

LEOData.ISS = LEOISSData;
LEOData.S = LEOSData;
LEOData.Other = LEOOtherData;



save('LEOData.mat', 'LEOData');    
save('totalLaunchData.mat','totalLaunchData');
save('LEOLaunchData.mat','LEOLaunchData');