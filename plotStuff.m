nep = readmatrix("neptune.csv");
plu = readmatrix("pluto.csv");
gm0 = 39.476926421372994; %1.327124400419393e11;
dt = 100000 / size(nep, 1);
ts = 0.0:dt:dt*(size(nep,1)-1);

nepElements = oak.tbp.elementsFromState(gm0, nep);
pluElements = oak.tbp.elementsFromState(gm0, plu);


figure(); hold on;
title("Neptune")
plot(ts, rad2deg(nepElements(:, 3) - nepElements(1, 3)))
xlabel("Time [years]")
ylabel("$\Delta i$ [deg]")

figure(); hold on;
title("Pluto")
plot(rad2deg(pluElements(:, 3) - pluElements(1, 3)))
xlabel("Time [years]")
ylabel("$\Delta i$ [deg]")

return
figure(); hold on;
% p1 = plot3([], [], [], LineWidth=3);
for k = 1:size(nepElements, 1)
    period = oak.tbp.period(gm0, pluElements(k, 1));
    ts = linspace(0, 2 * period, 200);
    nepStates = keplerPropagate(gm0, ts, nep(k, :));
    pluStates = keplerPropagate(gm0, ts, plu(k, :));
    rotMats = rotationMatrices(nepStates);

    pluRotated = ipermute(pagemtimes(pagetranspose(rotMats), permute(pluStates(:, 1:3), [2, 3, 1])), [2, 3, 1]);

%     set(p1, "XData", pluRotated(:, 1));
%     set(p1, "YData", pluRotated(:, 2));
%     set(p1, "ZData", pluRotated(:, 3));

clf

hold on;
%     plot3(nepStates(:, 1), nepStates(:, 2), nepStates(:,3 ));
%     plot3(pluStates(:, 1), pluStates(:, 2), pluStates(:,3 ));
    plot3(pluRotated(:, 1), pluRotated(:, 2), pluRotated(:,3 ));

    axis equal
    axis tight
    drawnow
end