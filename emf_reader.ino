#include <BasicLinearAlgebra.h>
#include <ElementStorage.h>
#include <Wire.h>
#include <Adafruit_Sensor.h>
#include <Adafruit_HMC5883_U.h>

using Vec3 = BLA::Matrix<3>;
using Mat3 = BLA::Matrix<3,3>;

#define ENABLE_TONE 1

#define DEBUG_PRINT 1
#define DEBUG_PRINT_CALIBRATION 0

#if DEBUG_PRINT
#define DEBUG_LOG(msg) Serial.println(msg)
#define DEBUG_LOGV(msg, i) do { Serial.print(msg); Serial.println(i) ; } while(0)
#else
#define DEBUG_LOG(msg)
#define DEBUG_LOGV(msg, i)
#endif

/* Assign a unique ID to this sensor at the same time */
Adafruit_HMC5883_Unified mag = Adafruit_HMC5883_Unified(12345);

void setup() {
  Serial.begin(9600);

  // Mag
  sensor_t sensor;
  mag.getSensor(&sensor);
  if (!mag.begin()) {
    while (true) {
      delay(1000);
    }
  }

#if DEBUG_PRINT_CALIBRATION
  Serial.println(F("Sensor Lab - IMU Calibration!"));
  Serial.println("Calibration filesys test");
  Serial.println("Looking for a magnetometer");
  Serial.println("Found addr 0x1c");
  Serial.println("Found a HMC5883 IMU");
  Serial.println("------------------------------------");
  mag.printSensorDetails();
  Serial.println("------------------------------------");
  Serial.println("");
  Serial.println("Looking for a gyroscope");
  Serial.println(F("Could not find a gyroscope, skipping!"));
  Serial.println("Looking for a accelerometer");
  Serial.println(F("Could not find a accelerometer, skipping!"));
#endif
  
  // Buzzer
  pinMode(11, OUTPUT);

  // LEDs
  pinMode(2, OUTPUT);
  pinMode(3, OUTPUT);
  pinMode(4, OUTPUT);
  pinMode(5, OUTPUT);
  pinMode(6, OUTPUT);
  digitalWrite(6, HIGH);

  // Button
  pinMode(7, OUTPUT);
  pinMode(8, INPUT);
  digitalWrite(7, HIGH);

  // "Boot" tone
#if ENABLE_TONE
  tone(11, 2000);
  delay(200);
  noTone(11);
#endif
}

void setEmf(int value)
{
  value = max(value, 1);
  value = min(value, 5);
  
  static int lastValue = 1;
  while (lastValue < value) {
    lastValue++;
    digitalWrite(7 - lastValue, HIGH);
  }
  while (lastValue > value) {
    digitalWrite(7 - lastValue, LOW);
    lastValue--;
  }

#if ENABLE_TONE
  if (value > 1) {
    tone(11, 696);
  } else {
    noTone(11);
  }
#endif
}

void modeDisplay(int bits)
{
  DEBUG_LOGV("display ", bits);
  for (int i = 0; i < 5; ++i)
  {
    digitalWrite(2 + i, bits & (1 << i));
  }
}

void debugHeading(float a, float b, float offPlane)
{
  float angle = atan2(a, b);
  int mid = offPlane >= 0.0f ? (1<<2) : 0;
  int mag = 1 + (int)(abs(angle) / (PI/3.0f));
  int magRight = ((mag & 1) << 1) | (mag >> 1);
  int magLeft = mag << 3;
  DEBUG_LOGV("mid", mid);
  DEBUG_LOGV("magRight", magRight);
  DEBUG_LOGV("magLeft", magLeft);
  modeDisplay(mid | (angle > 0.0f ? magRight : magLeft));
}

float randomf()
{
  static const long LONG_MAX = 2147483647;
  return random(LONG_MAX) / (float)LONG_MAX;
}

Vec3 randomSphere()
{
  float u = randomf();
  float v = randomf();
  float l = acos(2.0f * u - 1.0f) - PI / 2.0f;
  float t = 2.0f * PI * v;
  return Vec3(cos(l) * cos(t), cos(l) * sin(t), sin(l));
}

float length(Vec3 v)
{
  return sqrt(v(0)*v(0) + v(1)*v(1) + v(2)*v(2));
}

void normalize(Vec3& v)
{
  float l = length(v);
  v(0) /= l;
  v(1) /= l;
  v(2) /= l;
}

template<class T>
T interp(const T& a, const T& b, float x)
{
  return a + (b - a) * x;
}

float dot(const Vec3& a, const Vec3& b)
{
  return ((~a) * b)(0);
}

bool buttonDown()
{
  return digitalRead(8) == HIGH;
}

bool buttonPress()
{
  static bool lastState = false;
  bool state = buttonDown();
  bool pressed = lastState && !state;
  lastState = state;
  if (pressed) {
    DEBUG_LOG("Pressed the button");
  }
  return pressed;
}

Vec3 getMag()
{
  sensors_event_t event; 
  mag.getEvent(&event);

#if DEBUG_PRINT_CALIBRATION
  Serial.print("Raw:0,0,0,0,0,0,");
  Serial.print(event.magnetic.x); Serial.print(",");
  Serial.print(event.magnetic.y); Serial.print(",");
  Serial.print(event.magnetic.z);
  Serial.print("\n");
#endif

#if 0
  // https://github.com/nliaudat/magnetometer_calibration/blob/main/calibrate.py
  // pretty much didn't work at all
  static const float calibCenter[] = {
    32.04076239, 37.97431851, -65.1762087
  };
  static const float calibSphere[] = {
    2.31681277, 0.21903865, 0.300846858,
    0.21903865, 2.86506899, 0.15303164,
    0.30084685, 0.15303164, 2.47527054
  };
#endif

#if 1
  // Result of scipy optimizing without limitation
  // Likely involves some crazy skew to get it spherical
  static const Vec3 calibCenter{
    36.12895861,  41.88782217, -70.27233492
  };
  static const Mat3 calibSphere{
    1.87981651, -0.27689242,  0.66843622,
    0.52935566,  2.34396793,  0.74367147,
    -0.45195821, -0.21727938,  1.77663299
  };
#endif

  /*
  float cx = event.magnetic.x - calibCenter[0];
  float cy = event.magnetic.y - calibCenter[1];
  float cz = event.magnetic.z - calibCenter[2];
  result[0] = cx * calibSphere[0] + cy * calibSphere[1] + cz * calibSphere[2];
  result[1] = cx * calibSphere[3] + cy * calibSphere[4] + cz * calibSphere[5];
  result[2] = cx * calibSphere[6] + cy * calibSphere[7] + cz * calibSphere[8];
  */
  Vec3 result{event.magnetic.x, event.magnetic.y, event.magnetic.z};
  result = calibSphere * (result - calibCenter);
  
  //Serial.print(result[0]); Serial.print(" "); Serial.print(result[1]); Serial.print(" "); Serial.println(result[2]);

#if DEBUG_PRINT
  Serial << result;
  Serial.println("");
#endif
  return result;
}

enum State {
  WAITING,
  CHANGE_MODE,
  EVENT_FORCED,
  EVENT_MAGNITUDE,
  EVENT_DIRECTION,
};

enum Mode {
  MODE_NORMAL,
  MODE_SENSITIVITY_0,
  MODE_SENSITIVITY_1,
  MODE_SENSITIVITY_2,
  MODE_TEST_XY,
  MODE_TEST_YZ,
  MODE_TEST_XZ,
  MODE_TEST_MAGNITUDE,
  MODE_TEST_MAGNITUDE_EXP,

  MODE_COUNT,
};

static const unsigned long SELECT_MODE_TIME = 3000;

State state = WAITING;
Mode mode = MODE_NORMAL;
int ticker = 0;
unsigned long lastWaitTime = 0;
float sensitivity = 1.0f;
float targetMagnitude;
float targetMagnitudeInterp;
float targetMagnitudeBand = 50.0f;
Vec3 targetDirection;
Vec3 targetDirectionInterp;
float targetDirectionConeAngle = 1.0f;
int emfEventLevel = 1;
unsigned long emfEventDuration = 1000;

int randomEmfLevel()
{
  int r = random(100);
  if (r < 80) {
    Vec3 mag = getMag();
    float l = length(mag);
    if (l < 120.0f) return 2;
    if (l < 200.0f) return 3;
    if (l > 300.0f) return 4;
    return 5;
  } else {
    r = random(100);
    if (r < 50) return 2;
    if (r < 75) return 3;
    if (r < 90) return 4;
    return 5;
  }
}

unsigned long randomEventTime()
{
  return 10000 + random(40000 * sensitivity);
}

void enterState(int newState)
{
  // Old state cleanup
  DEBUG_LOGV("Leaving state ", state);
  switch (state) {
  case CHANGE_MODE:
    break;
  case EVENT_FORCED:
  case EVENT_MAGNITUDE:
  case EVENT_DIRECTION:
    modeDisplay(1<<4);
    setEmf(1);
    break;
  }

  DEBUG_LOGV("Entering state ", newState);
  switch (newState) {
  case WAITING:
    modeDisplay(1<<4);
    setEmf(1);
    mode = MODE_NORMAL;
    break;
  case CHANGE_MODE:
    mode = (Mode)((int)mode + 1);
    if (mode >= (int)MODE_COUNT) {
      mode = MODE_NORMAL;
    }
    DEBUG_LOGV("Setting mode ", mode);
    break;
  case EVENT_FORCED:
    {
      emfEventLevel = randomEmfLevel();
      emfEventDuration = randomEventTime();
      setEmf(emfEventLevel);
      break;
    }
  case EVENT_MAGNITUDE:
    {
      Vec3 mag = getMag();
      float r = randomf() * 4.0 - 1.0;
      float l = length(mag);
      float o = r * sensitivity * 0.1;
      targetMagnitude = l * (1.0 + o);
      targetMagnitudeInterp = targetMagnitude + o * (randomf() * 0.2f - 0.1f);
      targetMagnitudeBand = o * (randomf() * 0.5f + 0.5f);
      emfEventLevel = randomEmfLevel();
      emfEventDuration = randomEventTime() * 2;
      break;
    }
  case EVENT_DIRECTION:
    {
      Vec3 mag = getMag();
      // TODO: make these closer to the current direction based on sensitivity
      targetDirection = randomSphere();
      targetDirectionInterp = randomSphere();
      emfEventLevel = randomEmfLevel();
      emfEventDuration = randomEventTime();
      break;
    }
  }

  state = (State)newState;
  lastWaitTime = millis();
}

void loop() {
  if (buttonPress()){
    DEBUG_LOG("Changing mode");
    enterState(CHANGE_MODE);
  }

  unsigned long now = millis();

  if (mode != MODE_NORMAL) {
    if (now - lastWaitTime < SELECT_MODE_TIME) {
      switch (mode) {
      case MODE_NORMAL: modeDisplay(0x0); break;
      case MODE_SENSITIVITY_0: modeDisplay(0x1); break;
      case MODE_SENSITIVITY_1: modeDisplay(0x2); break;
      case MODE_SENSITIVITY_2: modeDisplay(0x3); break;
      case MODE_TEST_XY: modeDisplay(0x5); break;
      case MODE_TEST_YZ: modeDisplay(0x6); break;
      case MODE_TEST_XZ: modeDisplay(0x7); break;
      case MODE_TEST_MAGNITUDE: modeDisplay(0x8); break;
      case MODE_TEST_MAGNITUDE_EXP: modeDisplay(0x9); break;
      }
      return;
    }
    switch (mode) {
    case MODE_NORMAL: break;
    case MODE_SENSITIVITY_0: sensitivity = 0.3f; break;
    case MODE_SENSITIVITY_1: sensitivity = 1.0f; break;
    case MODE_SENSITIVITY_2: sensitivity = 3.0f; break;
    case MODE_TEST_XY:
    case MODE_TEST_YZ:
    case MODE_TEST_XZ:
      {
        Vec3 mag = getMag();
        switch (mode) {
          case MODE_TEST_XY: debugHeading(mag(0), mag(1), mag(2)); break;
          case MODE_TEST_YZ: debugHeading(mag(1), mag(2), mag(0)); break;
          case MODE_TEST_XZ: debugHeading(mag(0), mag(2), mag(1)); break;
        }
        // Keep this mode
        return;
      }
    case MODE_TEST_MAGNITUDE:
      {
        Vec3 mag = getMag();
        float l = length(mag);
        DEBUG_LOGV("Raw magnitude: ", l);
        modeDisplay(constrain((int)(sensitivity * l) / 10, 0, 31));
        l /= 4.0f;
        l *= l;
        tone(11, l);
        // Keep this mode
        return;
      }
    case MODE_TEST_MAGNITUDE_EXP:
      {
        Vec3 mag = getMag();
        float l = length(mag);
        DEBUG_LOGV("Raw magnitude: ", l);
        int v = 1<<4;
        if (l > 110.0f) v |= 1<<3;
        if (l > 150.0f) v |= 1<<2;
        if (l > 200.0f) v |= 1<<1;
        if (l > 300.0f) v |= 1<<0;
        modeDisplay(v);
        // Keep this mode
        return;
      }
    default:
      DEBUG_LOGV("Unhandled mode ", mode);
    }
    DEBUG_LOGV("Leaving mode select ", mode);
    enterState(WAITING);
  }
  
  switch (state) {
  case WAITING:
    if (now - lastWaitTime > 5000 / sensitivity) {
      int x = random(100);
      DEBUG_LOGV("Rolling 1d100: ", x);
      if (x < 2) {
        enterState(EVENT_FORCED);
      } else if (x < 7) {
        enterState(EVENT_MAGNITUDE);
      } else if (x < 10) {
        enterState(EVENT_DIRECTION);
      }
      lastWaitTime = now;
    }
    return;
  case EVENT_FORCED:
    if (now - lastWaitTime > emfEventDuration) {
      enterState(WAITING);
    }
    break;
  case EVENT_MAGNITUDE:
    {
      Vec3 mag = getMag();
      float l = length(mag);
      float t = interp(targetMagnitude, targetMagnitudeInterp, (now - lastWaitTime) / (float)emfEventDuration);
      DEBUG_LOGV("Compare magnitude ", t - l);
      if (l > t - targetMagnitudeBand && l < t + targetMagnitudeBand) {
        setEmf(emfEventLevel);
      } else {
        setEmf(1);
      }
      if (now - lastWaitTime > emfEventDuration) {
        enterState(WAITING);
      }
      break;
    }
  case EVENT_DIRECTION:
    {
      Vec3 mag = getMag();
      normalize(mag);
      Vec3 t = interp(targetDirection, targetDirectionInterp, (now - lastWaitTime) / (float)emfEventDuration);
      normalize(t);
      float a = acos(dot(mag, t));
      DEBUG_LOGV("Compare angle ", targetDirectionConeAngle - a);
      if (a < targetDirectionConeAngle) {
        setEmf(emfEventLevel);
      } else {
        setEmf(1);
      }
      if (now - lastWaitTime > emfEventDuration) {
        enterState(WAITING);
      }
      break;
    }
  }
}
