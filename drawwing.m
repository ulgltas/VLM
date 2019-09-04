clear all
close all

outfile
it=it+1;

if exist('vtail')
   vtail2.faces=vtail.faces(1:length(vtail.faces)/2,:);
   vtail2.vertices=vtail.vertices(1:length(vtail.vertices)/2,:);
   rudder2.faces=rudder.faces(1:length(rudder.faces)/2,:);
   rudder2.vertices=rudder.vertices(1:length(rudder.vertices)/2,:);
end

figure(1)
patch(wing,'FaceColor','b','EdgeColor','w');
patch(flap,'FaceColor','r','EdgeColor','w');
patch(aileron,'FaceColor','g','EdgeColor','w');
if exist('htail')
    patch(htail,'FaceColor','m','EdgeColor','w');
    patch(elevator,'FaceColor','c','EdgeColor','w');
end
if exist('vtail')
    patch(vtail2,'FaceColor','m','EdgeColor','w');
    patch(rudder2,'FaceColor','c','EdgeColor','w');
end
hold on
dummy=0;
for i=1:pwing.nwakes
    surf(pwing.xw(dummy+1:(dummy+pwing.wakelengths(i)+1),1:it+1), ...
        pwing.yw(dummy+1:(dummy+pwing.wakelengths(i)+1),1:it+1), ...
        pwing.zw(dummy+1:(dummy+pwing.wakelengths(i)+1),1:it+1),'FaceColor','none')
    dummy=dummy+pwing.wakelengths(i)+1;
end
dummy=0;
for i=1:pflap.nwakes
    surf(pflap.xw(dummy+1:(dummy+pflap.wakelengths(i)+1),1:it+1), ...
        pflap.yw(dummy+1:(dummy+pflap.wakelengths(i)+1),1:it+1), ...
        pflap.zw(dummy+1:(dummy+pflap.wakelengths(i)+1),1:it+1),'FaceColor','none')
    dummy=dummy+pflap.wakelengths(i)+1;
end
dummy=0;
for i=1:paileron.nwakes
    surf(paileron.xw(dummy+1:(dummy+paileron.wakelengths(i)+1),1:it+1), ...
        paileron.yw(dummy+1:(dummy+paileron.wakelengths(i)+1),1:it+1), ...
        paileron.zw(dummy+1:(dummy+paileron.wakelengths(i)+1),1:it+1),'FaceColor','none')
    dummy=dummy+paileron.wakelengths(i)+1;
end
if exist('htail')
    dummy=0;
    for i=1:phtail.nwakes
        surf(phtail.xw(dummy+1:(dummy+phtail.wakelengths(i)+1),1:it+1), ...
            phtail.yw(dummy+1:(dummy+phtail.wakelengths(i)+1),1:it+1), ...
            phtail.zw(dummy+1:(dummy+phtail.wakelengths(i)+1),1:it+1),'FaceColor','none')
        dummy=dummy+phtail.wakelengths(i)+1;
    end
    dummy=0;
    for i=1:pelevator.nwakes
        surf(pelevator.xw(dummy+1:(dummy+pelevator.wakelengths(i)+1),1:it+1), ...
            pelevator.yw(dummy+1:(dummy+pelevator.wakelengths(i)+1),1:it+1), ...
            pelevator.zw(dummy+1:(dummy+pelevator.wakelengths(i)+1),1:it+1),'FaceColor','none')
        dummy=dummy+pelevator.wakelengths(i)+1;
    end
end
if exist('vtail')
    dummy=0;
    for i=1:pvtail.nwakes-1
        surf(pvtail.xw(dummy+1:(dummy+pvtail.wakelengths(i)+1),1:it+1), ...
            pvtail.yw(dummy+1:(dummy+pvtail.wakelengths(i)+1),1:it+1), ...
            pvtail.zw(dummy+1:(dummy+pvtail.wakelengths(i)+1),1:it+1),'FaceColor','none')
        dummy=dummy+pvtail.wakelengths(i)+1;
    end
    dummy=0;
    for i=1:prudder.nwakes-1
        surf(prudder.xw(dummy+1:(dummy+prudder.wakelengths(i)+1),1:it+1), ...
            prudder.yw(dummy+1:(dummy+prudder.wakelengths(i)+1),1:it+1), ...
            prudder.zw(dummy+1:(dummy+prudder.wakelengths(i)+1),1:it+1),'FaceColor','none')
        dummy=dummy+prudder.wakelengths(i)+1;
    end
end
axis equal
xlabel('x')
ylabel('y')
zlabel('z')
view(3)

% figure(2)
% patch(wing,'FaceColor','none','EdgeColor','b');
% patch(flap,'FaceColor','none','EdgeColor','r');
% patch(aileron,'FaceColor','none','EdgeColor','g');
% if exist('htail')
%     patch(htail,'FaceColor','none','EdgeColor','m');
%     patch(elevator,'FaceColor','none','EdgeColor','c');
% end
% if exist('vtail')
%     patch(vtail,'FaceColor','none','EdgeColor','m');
%     patch(rudder,'FaceColor','none','EdgeColor','c');
% end
% axis equal
% xlabel('x')
% ylabel('m')
% zlabel('z')
% hold on
% plot3(wingcp(:,1),wingcp(:,2),wingcp(:,3),'bx','markersize',10)
% plot3(wingcv(:,1),wingcv(:,2),wingcv(:,3),'bo','markersize',10)
% plot3(flapcp(:,1),flapcp(:,2),flapcp(:,3),'rx','markersize',10)
% plot3(flapcv(:,1),flapcv(:,2),flapcv(:,3),'ro','markersize',10)
% plot3(aileroncp(:,1),aileroncp(:,2),aileroncp(:,3),'gx','markersize',10)
% plot3(aileroncv(:,1),aileroncv(:,2),aileroncv(:,3),'go','markersize',10)
% plot3(pwing.xw(:,1),pwing.yw(:,1),pwing.zw(:,1),'dk')
% plot3(pflap.xw(:,1),pflap.yw(:,1),pflap.zw(:,1),'^k')
% plot3(paileron.xw(:,1),paileron.yw(:,1),paileron.zw(:,1),'vk')
% if exist('htail')
%     plot3(htailcp(:,1),htailcp(:,2),htailcp(:,3),'mx','markersize',10)
%     plot3(htailcv(:,1),htailcv(:,2),htailcv(:,3),'mo','markersize',10)
%     plot3(phtail.xw(:,1),phtail.yw(:,1),phtail.zw(:,1),'*k')
%     plot3(elevatorcp(:,1),elevatorcp(:,2),elevatorcp(:,3),'cx','markersize',10)
%     plot3(elevatorcv(:,1),elevatorcv(:,2),elevatorcv(:,3),'co','markersize',10)
%     plot3(pelevator.xw(:,1),pelevator.yw(:,1),pelevator.zw(:,1),'*k')
% end
% if exist('vtail')
%     plot3(vtailcp(:,1),vtailcp(:,2),vtailcp(:,3),'mx','markersize',10)
%     plot3(vtailcv(:,1),vtailcv(:,2),vtailcv(:,3),'mo','markersize',10)
%     plot3(pvtail.xw(:,1),pvtail.yw(:,1),pvtail.zw(:,1),'*k')
%     plot3(ruddercp(:,1),ruddercp(:,2),ruddercp(:,3),'cx','markersize',10)
%     plot3(ruddercv(:,1),ruddercv(:,2),ruddercv(:,3),'co','markersize',10)
%     plot3(prudder.xw(:,1),prudder.yw(:,1),prudder.zw(:,1),'*k')
% end
% 
% figure(3)
% patch(wing,'FaceColor','none','EdgeColor','b');
% patch(flap,'FaceColor','none','EdgeColor','r');
% patch(aileron,'FaceColor','none','EdgeColor','g');
% if exist('htail')
%     patch(htail,'FaceColor','none','EdgeColor','m');
%     patch(elevator,'FaceColor','none','EdgeColor','c');
% end
% if exist('vtail')
%     patch(vtail,'FaceColor','none','EdgeColor','m');
%     patch(rudder,'FaceColor','none','EdgeColor','c');
% end
% axis equal
% xlabel('x')
% ylabel('m')
% zlabel('z')
% hold on
% quiver3(wingcp(:,1),wingcp(:,2),wingcp(:,3),wingnorm(:,1),wingnorm(:,2),wingnorm(:,3),0)
% quiver3(aileroncp(:,1),aileroncp(:,2),aileroncp(:,3),aileronnorm(:,1),aileronnorm(:,2),aileronnorm(:,3),0)
% quiver3(flapcp(:,1),flapcp(:,2),flapcp(:,3),flapnorm(:,1),flapnorm(:,2),flapnorm(:,3),0)
% if exist('htail')
%     quiver3(elevatorcp(:,1),elevatorcp(:,2),elevatorcp(:,3),elevatornorm(:,1),elevatornorm(:,2),elevatornorm(:,3),0)
%     quiver3(htailcp(:,1),htailcp(:,2),htailcp(:,3),htailnorm(:,1),htailnorm(:,2),htailnorm(:,3),0)
% end
% if exist('vtail')
%     quiver3(ruddercp(:,1),ruddercp(:,2),ruddercp(:,3),ruddernorm(:,1),ruddernorm(:,2),ruddernorm(:,3),0)
%     quiver3(vtailcp(:,1),vtailcp(:,2),vtailcp(:,3),vtailnorm(:,1),vtailnorm(:,2),vtailnorm(:,3),0)
% end
% 
% figure(4)
% patch(wing,'FaceColor','none','EdgeColor','b');
% patch(flap,'FaceColor','none','EdgeColor','r');
% patch(aileron,'FaceColor','none','EdgeColor','g');
% if exist('htail')
%     patch(htail,'FaceColor','none','EdgeColor','m');
%     patch(elevator,'FaceColor','none','EdgeColor','c');
% end
% if exist('vtail')
%     patch(vtail,'FaceColor','none','EdgeColor','m');
%     patch(rudder,'FaceColor','none','EdgeColor','c');
% end
% axis equal
% xlabel('x')
% ylabel('m')
% zlabel('z')
% hold on
% quiver3(wingcp(:,1),wingcp(:,2),wingcp(:,3),wingtangx(:,1),wingtangx(:,2),wingtangx(:,3),0)
% quiver3(aileroncp(:,1),aileroncp(:,2),aileroncp(:,3),ailerontangx(:,1),ailerontangx(:,2),ailerontangx(:,3),0)
% quiver3(flapcp(:,1),flapcp(:,2),flapcp(:,3),flaptangx(:,1),flaptangx(:,2),flaptangx(:,3),0)
% if exist('htail')
%     quiver3(elevatorcp(:,1),elevatorcp(:,2),elevatorcp(:,3),elevatortangx(:,1),elevatortangx(:,2),elevatortangx(:,3),0)
%     quiver3(htailcp(:,1),htailcp(:,2),htailcp(:,3),htailtangx(:,1),htailtangx(:,2),htailtangx(:,3),0)
% end
% if exist('vtail')
%     quiver3(ruddercp(:,1),ruddercp(:,2),ruddercp(:,3),ruddertangx(:,1),ruddertangx(:,2),ruddertangx(:,3),0)
%     quiver3(vtailcp(:,1),vtailcp(:,2),vtailcp(:,3),vtailtangx(:,1),vtailtangx(:,2),vtailtangx(:,3),0)
% end
% 
% figure(5)
% patch(wing,'FaceColor','none','EdgeColor','b');
% patch(flap,'FaceColor','none','EdgeColor','r');
% patch(aileron,'FaceColor','none','EdgeColor','g');
% if exist('htail')
%     patch(htail,'FaceColor','none','EdgeColor','m');
%     patch(elevator,'FaceColor','none','EdgeColor','c');
% end
% if exist('vtail')
%     patch(vtail,'FaceColor','none','EdgeColor','m');
%     patch(rudder,'FaceColor','none','EdgeColor','c');
% end
% axis equal
% xlabel('x')
% ylabel('m')
% zlabel('z')
% hold on
% quiver3(wingcp(:,1),wingcp(:,2),wingcp(:,3),wingtangy(:,1),wingtangy(:,2),wingtangy(:,3),0)
% quiver3(aileroncp(:,1),aileroncp(:,2),aileroncp(:,3),ailerontangy(:,1),ailerontangy(:,2),ailerontangy(:,3),0)
% quiver3(flapcp(:,1),flapcp(:,2),flapcp(:,3),flaptangy(:,1),flaptangy(:,2),flaptangy(:,3),0)
% if exist('htail')
%     quiver3(elevatorcp(:,1),elevatorcp(:,2),elevatorcp(:,3),elevatortangy(:,1),elevatortangy(:,2),elevatortangy(:,3),0)
%     quiver3(htailcp(:,1),htailcp(:,2),htailcp(:,3),htailtangy(:,1),htailtangy(:,2),htailtangy(:,3),0)
% end
% if exist('vtail')
%     quiver3(ruddercp(:,1),ruddercp(:,2),ruddercp(:,3),ruddertangy(:,1),ruddertangy(:,2),ruddertangy(:,3),0)
%     quiver3(vtailcp(:,1),vtailcp(:,2),vtailcp(:,3),vtailtangy(:,1),vtailtangy(:,2),vtailtangy(:,3),0)
% end
% 
% figure(6)
% patch(wing,'FaceColor','none','EdgeColor','b');
% patch(flap,'FaceColor','none','EdgeColor','r');
% patch(aileron,'FaceColor','none','EdgeColor','g');
% if exist('htail')
%     patch(htail,'FaceColor','none','EdgeColor','m');
%     patch(elevator,'FaceColor','none','EdgeColor','c');
% end
% if exist('vtail')
%     patch(vtail,'FaceColor','none','EdgeColor','m');
%     patch(rudder,'FaceColor','none','EdgeColor','c');
% end
% axis equal
% xlabel('x')
% ylabel('m')
% zlabel('z')
% hold on
% quiver3(wingcp(:,1),wingcp(:,2),wingcp(:,3),pwing.uvw(:,1),pwing.uvw(:,2),pwing.uvw(:,3),0)
% quiver3(aileroncp(:,1),aileroncp(:,2),aileroncp(:,3),paileron.uvw(:,1),paileron.uvw(:,2),paileron.uvw(:,3),0)
% quiver3(flapcp(:,1),flapcp(:,2),flapcp(:,3),pflap.uvw(:,1),pflap.uvw(:,2),pflap.uvw(:,3),0)
% if exist('htail')
%     quiver3(htailcp(:,1),htailcp(:,2),htailcp(:,3),phtail.uvw(:,1),phtail.uvw(:,2),phtail.uvw(:,3),0)
%     quiver3(elevatorcp(:,1),elevatorcp(:,2),elevatorcp(:,3),pelevator.uvw(:,1),pelevator.uvw(:,2),pelevator.uvw(:,3),0)
% end
% if exist('vtail')
%     quiver3(vtailcp(:,1),vtailcp(:,2),vtailcp(:,3),pvtail.uvw(:,1),pvtail.uvw(:,2),pvtail.uvw(:,3),0)
%     quiver3(ruddercp(:,1),ruddercp(:,2),ruddercp(:,3),prudder.uvw(:,1),prudder.uvw(:,2),prudder.uvw(:,3),0)
% end

figure(8)
plot((0:it-1)*dt,totalforce([4 2 3],1:it))
xlabel('time (s)')
ylabel('force (N)')
legend('f_x','f_y','f_z')

% Corrseponds to points wingcp
Wing_force_x=wingDeltap.*wingnsurf.*wingnorm(:,1);
Wing_force_y=wingDeltap.*wingnsurf.*wingnorm(:,2);
Wing_force_z=wingDeltap.*wingnsurf.*wingnorm(:,3);

figure(103)
patch(wing,'FaceColor','none','EdgeColor','b');
patch(flap,'FaceColor','none','EdgeColor','r');
patch(aileron,'FaceColor','none','EdgeColor','g');
if exist('htail')
    patch(htail,'FaceColor','none','EdgeColor','m');
    patch(elevator,'FaceColor','none','EdgeColor','c');
end
if exist('vtail')
    patch(vtail2,'FaceColor','none','EdgeColor','m');
    patch(rudder2,'FaceColor','none','EdgeColor','c');
end
axis equal
xlabel('x')
ylabel('m')
zlabel('z')
hold on
scalefact=10;
wingforces=wingDeltap.*wingnsurf/scalefact;
aileronforces=aileronDeltap.*aileronnsurf/scalefact;
flapforces=flapDeltap.*flapnsurf/scalefact;
htailforces=htailDeltap.*htailnsurf/scalefact;
elevatorforces=elevatorDeltap.*elevatornsurf/scalefact;
vtailforces=vtailDeltap(1:length(vtail.faces)/2).*vtailnsurf(1:length(vtail.faces)/2)/scalefact;
rudderforces=rudderDeltap(1:length(rudder.faces)/2).*ruddernsurf(1:length(rudder.faces)/2)/scalefact;

quiver3(wingcp(:,1),wingcp(:,2),wingcp(:,3),wingnorm(:,1).*wingforces,wingnorm(:,2).*wingforces,wingnorm(:,3).*wingforces,0)
if ~isempty(aileron.faces)
    quiver3(aileroncp(:,1),aileroncp(:,2),aileroncp(:,3),aileronnorm(:,1).*aileronforces,aileronnorm(:,2).*aileronforces,aileronnorm(:,3).*aileronforces,0)
end
if ~isempty(flap.faces)
    quiver3(flapcp(:,1),flapcp(:,2),flapcp(:,3),flapnorm(:,1).*flapforces,flapnorm(:,2).*flapforces,flapnorm(:,3).*flapforces,0)
end
if ~isempty(elevator.faces)
    quiver3(elevatorcp(:,1),elevatorcp(:,2),elevatorcp(:,3),elevatornorm(:,1).*elevatorforces,elevatornorm(:,2).*elevatorforces,elevatornorm(:,3).*elevatorforces,0)
end
if ~isempty(htail.faces)
    quiver3(htailcp(:,1),htailcp(:,2),htailcp(:,3),htailnorm(:,1).*htailforces,htailnorm(:,2).*htailforces,htailnorm(:,3).*htailforces,0)
end
if ~isempty(rudder.faces)
    quiver3(ruddercp((1:length(rudder.faces)/2),1),ruddercp((1:length(rudder.faces)/2),2),ruddercp((1:length(rudder.faces)/2),3),ruddernorm((1:length(rudder.faces)/2),1).*rudderforces,ruddernorm((1:length(rudder.faces)/2),2).*rudderforces,ruddernorm((1:length(rudder.faces)/2),3).*rudderforces,0)
end
if ~isempty(vtail.faces)
    quiver3(vtailcp((1:length(vtail.faces)/2),1),vtailcp((1:length(vtail.faces)/2),2),vtailcp((1:length(vtail.faces)/2),3),vtailnorm((1:length(vtail.faces)/2),1).*vtailforces,vtailnorm((1:length(vtail.faces)/2),2).*vtailforces,vtailnorm((1:length(vtail.faces)/2),3).*vtailforces,0)
end
title('Aerodynamic forces normal to each panel')
