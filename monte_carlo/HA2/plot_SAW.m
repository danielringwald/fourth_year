nodes = [0 0; 0 1; 1 1; 2 1; 2 0; 2 -1; 1 -1; 1 0];
x1 = -1; y1 = -2;
x2 = 3;  y2 = 3;
figure
hold on
plot(nodes(:,1),nodes(:,2),'--*')
axis([x1 x2 y1 y2])
title("Spiraling walk")
xlabel("x")
ylabel("y")
grid on
hold off
