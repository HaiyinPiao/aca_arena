/** Example 004 Movement

This Tutorial shows how to move and animate SceneNodes. The
basic concept of SceneNodeAnimators is shown as well as manual
movement of nodes using the keyboard.  We'll demonstrate framerate
independent movement, which means moving by an amount dependent
on the duration of the last run of the Irrlicht loop.

Example 19.MouseAndJoystick shows how to handle those kinds of input.

As always, I include the header files, use the irr namespace,
and tell the linker to link with the .lib file.
*/
#ifdef _MSC_VER
// We'll also define this to stop MSVC complaining about sprintf().
#define _CRT_SECURE_NO_WARNINGS
#pragma comment(lib, "Irrlicht.lib")
#endif

#include <irrlicht.h>
#include "driverChoice.h"
#include <cstdio>
#include <ctime>
#include <iostream>

extern "C"
{
#include "f16-hrt/f16.h"
}

using namespace std;
using namespace irr;

// integration steps, in ms.
const float g_intstep=30;

float g_elpstime=0;

class MyEventReceiver : public IEventReceiver
{
public:
	// This is the one method that we have to implement
	virtual bool OnEvent(const SEvent& event)
	{
		// Remember whether each key is down or up
		if (event.EventType == irr::EET_KEY_INPUT_EVENT)
			KeyIsDown[event.KeyInput.Key] = event.KeyInput.PressedDown;

		return false;
	}

	// This is used to check whether a key is being held down
	virtual bool IsKeyDown(EKEY_CODE keyCode) const
	{
		return KeyIsDown[keyCode];
	}
	
	MyEventReceiver()
	{
		for (u32 i=0; i<KEY_KEY_CODES_COUNT; ++i)
			KeyIsDown[i] = false;
	}

private:
	// We use this array to store the current state of each key
	bool KeyIsDown[KEY_KEY_CODES_COUNT];
};

int main()
{
    f16_initialize();

    // input
    f16_U.n_nc = 0.0f;
    f16_U.n_xc = 0.0f;
    f16_U.mu_dotc = 3.14f/2.0f;

    // holds the simulation steps
    unsigned int k=0;

	MyEventReceiver receiver;

    IrrlichtDevice* device = createDevice(video::EDT_OPENGL,
            core::dimension2d<u32>(1440, 960), 16, false, false, false, &receiver);

	if (device == 0)
		return 1; // could not create selected driver.

	video::IVideoDriver* driver = device->getVideoDriver();
	scene::ISceneManager* smgr = device->getSceneManager();

	scene::IAnimatedMeshSceneNode* anms =
        smgr->addAnimatedMeshSceneNode(smgr->getMesh("../aca_arena/T50-airframe.3ds"));

	if (anms)
	{
        anms->setMaterialTexture(0, driver->getTexture("../aca_arena/T50-airframe.png"));
        /*scene::ISceneNodeAnimator* anim =
            smgr->createFlyCircleAnimator(core::vector3df(0,0,1), 10.0f);
        if (anim)
        {
            anms->addAnimator(anim);
            anim->drop();
        }*/

        anms->setMaterialFlag(video::EMF_LIGHTING, false);

		anms->setFrameLoop(0, 13);
		anms->setAnimationSpeed(15);

        anms->setScale(core::vector3df(300.f,300.f,300.f));
        /*anms->setRotation(core::vector3df(0,0,0));
        anms->setPosition(core::vector3df(0,0,0));*/
//		anms->setMaterialTexture(0, driver->getTexture("../../media/sydney.bmp"));

	}

	/*
	To be able to look at and move around in this scene, we create a first
	person shooter style camera and make the mouse cursor invisible.
	*/
    //smgr->addCameraSceneNodeFPS();
    smgr->addCameraSceneNode(0, core::vector3df(0,-1500,-1000), core::vector3df(0,1000,1000));
    //device->getCursorControl()->setVisible(false);

	int lastFPS = -1;

	// In order to do framerate independent movement, we have to know
	// how long it was since the last frame
	u32 then = device->getTimer()->getTime();

	// This is the movemen speed in units per second.
	const f32 MOVEMENT_SPEED = 5.f;

	while(device->run())
	{
		// Work out a frame delta time.
		const u32 now = device->getTimer()->getTime();
		const f32 frameDeltaTime = (f32)(now - then) / 1000.f; // Time in seconds
		then = now;

        g_elpstime += frameDeltaTime;

        if( g_elpstime>(g_intstep/1000.0f) )
        {
            //output
            float x=f16_Y.Xg[0];
            float y=f16_Y.Xg[1];
            float z=f16_Y.Xg[2];

            float mu=f16_Y.attitude_g[0];
            float gamma=f16_Y.attitude_g[1];
            float psi=f16_Y.attitude_g[2];

            f16_step();

            anms->setRotation(core::vector3df(mu*57.3+180,gamma*57.3,psi*57.3));
            anms->setPosition(core::vector3df(x,y,z+1000));

            g_elpstime=0;
        }

        driver->beginScene(true, true, video::SColor(255,255,255,255));

		smgr->drawAll(); // draw the 3d scene
		device->getGUIEnvironment()->drawAll(); // draw the gui environment (the logo)

		driver->endScene();

		int fps = driver->getFPS();

		if (lastFPS != fps)
		{
            core::stringw tmp(L"ACA::Arena - haiyinpiao@qq.com [");
			tmp += driver->getName();
			tmp += L"] fps: ";
			tmp += fps;

			device->setWindowCaption(tmp.c_str());
			lastFPS = fps;
		}
	}

	/*
	In the end, delete the Irrlicht device.
	*/
	device->drop();

    f16_terminate();
	
	return 0;
}

/*
That's it. Compile and play around with the program.
**/
