import asyncio
import struct

async def connect():
    reader, writer = await asyncio.open_connection('::1', 1383)

    writer.write(b'ARTIQ moninj\n'.encode())
    endian_tmp = await self.reader.read(1)
    endian = None
    if endian_tmp == b'e':
        endian = '<'
    elif endian_tmp == b'E':
        endian = '>'
    else:
        print('ERROR')
    await writer.drain()
    print(endian)
    # packet = struct.pack(endian + "bblb", 0, True, channel, probe)
    # writer.write(packet)
    #
    # data = await reader.read(100)
    # print(f'Received: {data.decode()!r}')

    print('Close the connection')
    writer.close()
    await writer.wait_closed()

asyncio.run(connect())
