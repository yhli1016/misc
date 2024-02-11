#! /usr/bin/env python


def print_banner_line(text: str,
                      width: int = 80,
                      mark: str = "-",
                      end: str = "#") -> None:
    """
    Print a banner like '#--------------- FOO ---------------#' to stdout.

    :param text: central text in the banner
    :param width: total width of the banner
    :param mark: marker character
    :param end: end character prepended and appended to the banner
    :return: None
    """
    num_marks_total = width - len(text) - 4
    num_marks_left = num_marks_total // 2
    num_marks_right = num_marks_total - num_marks_left
    banner_with_marks = end + mark * num_marks_left
    banner_with_marks += " %s " % text
    banner_with_marks += mark * num_marks_right + end
    print(banner_with_marks)


def print_banner_block(text: str,
                       width: int = 80,
                       mark: str = "-",
                       end: str = "#") -> None:
    """
    Print a banner like
    #----------------------------------#
    #               FOO                #
    #----------------------------------#
    to stdout.

    :param text: string, central text in the banner
    :param width: integer, total width of the banner
    :param mark: character, marker
    :param end: character, end prepended and appended to the banner
    :return: None
    """
    num_spaces_total = width - len(text) - 2
    num_spaces_left = num_spaces_total // 2
    num_spaces_right = num_spaces_total - num_spaces_left
    banner_with_spaces = end + " " * num_spaces_left
    banner_with_spaces += text
    banner_with_spaces += " " * num_spaces_right + end
    border = end + mark * (width - 2) + end
    print(border)
    print(banner_with_spaces)
    print(border)


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("text", type=str,  action="store")
    parser.add_argument("--width", "-w",   type=int, action="store",
                        default="80")
    parser.add_argument("--type", "-t", type=str, action="store",
                        default="line")
    parser.add_argument("--mark", "-m", type=str, action="store",
                        default="-")
    parser.add_argument("--end", "-e", type=str, action="store",
                        default="#")
    args = parser.parse_args()
    if args.type == "line":
        print_banner_line(args.text, args.width, args.mark, args.end)
    elif args.type == "block":
        print_banner_block(args.text, args.width, args.mark, args.end)
    else:
        raise ValueError("ERROR: Unknown banner type '%s'" % args.flavor)
