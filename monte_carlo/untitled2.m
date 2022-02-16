x1 = -10; y1 = -10;
x2 = 10;  y2 = 10;
figure
hold on
plot(X(:,1,1),X(:,2,1),'--*')
axis([x1 x2 y1 y2])
grid on
hold off
figure
hold on
plot(X(:,1,2),X(:,2,2),'--*')
axis([x1 x2 y1 y2])
grid on
hold off
figure
hold on
plot(X(:,1,3),X(:,2,3),'--*')
axis([x1 x2 y1 y2])
grid on
hold off

for i=4:6
    figure
    hold on
    plot(X(:,1,i),X(:,2,i),'--*')
    axis([x1 x2 y1 y2])
    grid on
    hold off

end