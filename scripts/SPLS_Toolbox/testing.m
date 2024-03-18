%% testing

pd = fitdist(RHO_b_collection,'Normal')
area(RHO_b_collection)
diff(pd([-2.5 2.5]))
plot(RHO_b_collection)
trapz(RHO_b_collection)
trapz(0.25, RHO_b_collection)

X = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85];
rA = [-.0053, -.0052, -.0050, -.0045, -.0040, -.0033, -.0025, -.0018, -.00125, -.0010];
trapz(X,-1./rA)
x=1:5000;
y=RHO_b_collection;
h=histogram(RHO_b_collection, 'BinWidth', 0.001)
trapz(h)
y = sort(RHO_b_collection);
x=1:5000;
unique(RHO_b_collection)

>> x=1:10
x =
     1     2     3     4     5     6     7     8     9    10
>> y=rand(1,10)
y =
  Columns 1 through 7
    0.8147    0.9058    0.1270    0.9134    0.6324    0.0975    0.2785
  Columns 8 through 10
    0.5469    0.9575    0.9649
>> trapz(x(4:6),y(4:6))
ans =
    1.1378
    
    X = 0:pi/100:pi;
Calculate the sine of X.

Y = sin(X);
Integrate Y using trapz.

Q = trapz(X(98:end),Y(98:end))

plot(X,Y)

bins = 0:0.001:1;

RHO_test = 0.2163;

fun = @(x_var, y_var)abs(x_var-y_var);
fun(RHO_test,3)
for i=1:size(bins,2)
   evalbins(i) = fun(bins(i), RHO_test); 
end

h = histogram(RHO_b_collection, 'BinWidth', 0.001);
x = 0:0.001:((h.NumBins-1)*0.001);
y = h.Values;

[~, index] = min(abs(x-RHO_test));
auc_ratio = trapz(x(index:end), y(index:end))/trapz(x,y)
auc_ratio = sum(cumtrapz(x(index:end), y(index:end)))/sum(cumtrapz(x,y))

IN.dist = RHO_b_collection;
IN.val = RHO_test;

[p_val, h] = dp_auc_testing(IN)
testing_nr=2

switch testing_nr
    case 2
        disp('1');
    case 1
        disp('2');
        
    case 3
        disp('3');
        
end
