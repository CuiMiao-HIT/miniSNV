from collections import namedtuple

CommandOption = namedtuple('CommandOption', ['option', 'value'])
CommandOptionWithNoValue = namedtuple('CommandOptionWithNoValue', ['option'])
ExecuteCommand = namedtuple('ExecuteCommand', ['bin', 'bin_value'])


def command_option_string_from(command):
    if isinstance(command, CommandOption):
        return "--{} {}".format(command.option, command.value) if command.value is not None else None
    elif isinstance(command, CommandOptionWithNoValue):
        return "--{}".format(command.option)
    elif isinstance(command, ExecuteCommand):
        return " ".join([command.bin, command.bin_value])
    else:
        return command


def command_string_from(command_options):
    return " ".join(x for x in map(command_option_string_from, command_options) if x is not None)

