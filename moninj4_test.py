"""
moninj: try directly talking to core
"""
import asyncio

from artiq.coredevice.comm_moninj import CommMonInj, TTLOverride


def monitor_cb(self, channel, probe, value):
    pass

def injection_status_cb(self, channel, override, value):
    pass

def main():
    loop = asyncio.get_event_loop()
    th1 = CommMonInj(monitor_cb, injection_status_cb)
    loop.run_until_complete(th1.connect('192.168.1.75', 1383))
    print(th1)
    th1.inject(8, TTLOverride.level.value, 1)
    th1.inject(8, TTLOverride.oe.value, 1)
    th1.inject(8, TTLOverride.en.value, 1)

if __name__ == "__main__":
    main()